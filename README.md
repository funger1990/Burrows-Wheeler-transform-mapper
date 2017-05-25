# Burrows-Wheeler-transform-mapper
super-fast search for sequence by Burrows-Wheeler-transform and FM-index


### required python modules
1. Numpy
2. Pandas
3. pyTable
4. Biopython

### reference data
please download [UniProt/Swiss-Prot](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) fasta.gz file

### note
tools_karkkainen_sanders.py is copied from [pysuffix](https://code.google.com/archive/p/pysuffix/) to build suffix array efficiently in O(n * log(n))
