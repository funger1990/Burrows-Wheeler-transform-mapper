#!/usr/bin/env python
from __future__ import division

__author__ = 'Fan Zhang'


"""
implement Burrows-Wheeler Transform
"""

from collections import Counter
import bisect
import sys
import gzip

import numpy as np
import pandas as pd
from Bio import SeqIO


# pysuffix
import tools_karkkainen_sanders as tks


# ----------------------------------------
def build_suffix_array_naive(t):
    # make suffix array for simplicity
    sa = sorted(range(len(t)), key=lambda x: t[x:])
    return sa

def build_suffix_array_pysuffix(t):
    sa =  tks.simple_kark_sort(t)[:len(t)]
    return sa


class BwtProtein(object):
    def __init__(self):
        self.alphabet = list('ACDEFGHIKLMNPQRSTVWY')

    def concat_fasta(self, fasta):
        # concatenate fasta to a long string, and generate index
        separator = '|'
        t = []
        pos = 0
        dict_protein = {}
        coord_protein = {}
        with gzip.open(fasta) as f:
            for record in SeqIO.parse(f, 'fasta'):
                coord_protein[pos] = record.id
                t.append(str(record.seq))
                pos += len(record.seq) + 1
                dict_protein[record.id] = str(record.seq)
        t = separator.join(t) + '$'
        coord_protein = pd.Series(coord_protein)
        return t, coord_protein, dict_protein

    def bwt_via_sa(self, t, func_build_sa=build_suffix_array_pysuffix):
        bw = []
        sa = func_build_sa(t)
        for i in sa:
            bw.append(t[i-1])
        bw = ''.join(bw)
        return bw, sa

    def create_first(self, bw):
        # make concise representation of first BWM column
        first = {}
        tot = 0
        for i, cnt in sorted(Counter(bw).items()):
            first[i] = tot
            tot += cnt
        return first

    def create_brank_checkpoint(self, bw, brank_bin):
        # calculate tally checkpoint of B-rank
        # count character before (exclusive)
        n = len(bw)
        counter_to_df = lambda x: pd.DataFrame.from_dict(x, orient='index').T
        counter = Counter({bw[0]: 0})
        df_list = [counter_to_df(counter)]

        for i in xrange(0, n, brank_bin):
            counter += Counter(bw[i:i+brank_bin])
            df_list.append(counter_to_df(counter))

        # by default sorted
        df = pd.concat(df_list)
        # bug!! must record '&'
        # df.drop('$', axis=1, inplace=True)
        df.fillna(0, inplace=True)
        # add total counts
        df.index = range(0, n, brank_bin) + [n]
        df = df.astype(np.uint32)

        return df

    def create_sa_checkpoint(self, sa, sa_bin):
        n = len(sa)
        sa_checkpoint = {}
        for i, x in enumerate(sa):
            if x % sa_bin == 0:
                sa_checkpoint[i] = x
        return sa_checkpoint

    def build_index(self, fasta, index, brank_bin=128, sa_bin=32):
        print 'concatenate fasta'
        t, coord_protein, dict_protein = self.concat_fasta(fasta)
        print 'Burrows-Wheeler transform via suffix array...'
        bw, sa = self.bwt_via_sa(t)
        print 'record checkpoint of B-rank and F column...'
        first = self.create_first(bw)
        brank_checkpoint = self.create_brank_checkpoint(bw, brank_bin)
        sa_checkpoint = self.create_sa_checkpoint(sa, sa_bin)
        data = pd.Series({'dict_protein': dict_protein, 'bw': bw, 'first': first, 'n': len(bw),
                          'brank_bin': brank_bin, 'sa_bin': sa_bin, 'sa_checkpoint': sa_checkpoint})

        # save as hdf
        print 'save index'
        with pd.HDFStore(index) as store:
            store['coord_protein'] = coord_protein
            store['data'] = data
            store['brank_checkpoint'] = brank_checkpoint


    def load_index(self, index):
        with pd.HDFStore(index) as store:
            self.coord_protein = store['coord_protein']
            self.brank_checkpoint = store['brank_checkpoint']
            self.dict_protein = store['data']['dict_protein']
            self.bw = store['data']['bw']
            self.n = store['data']['n']
            self.first = store['data']['first']
            self.brank_bin = store['data']['brank_bin']
            self.sa_bin = store['data']['sa_bin']
            self.sa_checkpoint = store['data']['sa_checkpoint']

    def lf_mapping(self, idx, c=None):
        # exact rather than nearest
        if c is None:
            c = self.bw[idx]
        # calculate B-rank with tally checkpoint
        brank_idx = idx // self.brank_bin * self.brank_bin
        brank = self.brank_checkpoint.loc[brank_idx, c] + self.bw[brank_idx : idx].count(c)
        # LF mapping
        idx = self.first[c] + brank
        return idx

    def query_bwt(self, pep):
        # right exclusive
        idx_lo = 0
        idx_hi = self.n
        # search from end to begin, LF mapping
        for c in pep[::-1]:
            idx_lo = self.lf_mapping(idx_lo, c)
            idx_hi = self.lf_mapping(idx_hi, c)
            if idx_lo == idx_hi:
                break

        # lookup suffix array checkpoint
        loc_list = []
        for idx in xrange(idx_lo, idx_hi):
            loc = 0
            while idx not in self.sa_checkpoint:
                loc += 1
                idx = self.lf_mapping(idx)
            loc += self.sa_checkpoint[idx]
            loc_list.append(loc)

        return loc_list

    def map_peptide(self, pep):
        # check pep
        pep = pep.upper()
        if set(list(pep)) - set(self.alphabet):
            print 'wrong peptide', pep
            return 0, '', '', 0, ''

        # map peptide to concatenated t via bwt
        loc_list = self.query_bwt(pep)

        # map to protein and find context
        protein_list = []
        pos_list = []
        context_list = []

        for loc in loc_list:
            protein_beg = self.coord_protein.index[bisect.bisect(self.coord_protein.index, loc) - 1]
            protein_name = self.coord_protein[protein_beg]
            pep_pos = loc - protein_beg
            upstream = self.dict_protein[protein_name][max(0, pep_pos-3) : pep_pos]
            downstream = self.dict_protein[protein_name][pep_pos + len(pep) : pep_pos + len(pep) + 3]
            pep_context = upstream + ':' + downstream
            protein_list.append(protein_name)
            pos_list.append(pep_pos)
            context_list.append(pep_context)

        n_protein = len(protein_list)
        protein_string = ','.join(protein_list)
        pos_string = ','.join(map(str, pos_list))
        context_list = list(set(context_list))
        n_context = len(context_list)
        context_string = ','.join(context_list)

        return n_protein, protein_string, pos_string, n_context, context_string


# ----------------------------------------
if __name__ == '__main__':
    print sys.platform

    if sys.platform == 'win32':
        fasta = 'G:/data/uniprot/uniprot_sprot_human_rename.fasta.gz'
        index = 'G:/data/uniprot/uniprot_bwt_index.h5'
    else:
        fasta = '/mnt/g/data/uniprot/uniprot_sprot_human_rename.fasta.gz'
        index = '/mnt/g/data/uniprot/uniprot_bwt_index.h5'

    B = BwtProtein()
    # B.build_index(fasta, index)
    B.load_index(index)
    print B.map_peptide('AERYDDMAA')

