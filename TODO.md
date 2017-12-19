# list of TODOs for sesbio
## gene_annotation
script: get_intervals.pl
 - [ ] find replacement for Set::IntervalTree, which is not compiling (correction: with some c++ compilers).
script: seqio_benchmarks.pl
 - [ ] add transposome seqfactory benchmarks, which are 3X faster than the transposome seqio methods
script: clean_multifasta.pl
 - [ ] make script work with protein sequences

## genome_assembly
script: gsAssembler.pl
 - [ ] remove hard coded paths and implement some exception handling for command execution

## transposon_annotation
script: find_graph_edges.pl
 - [ ] remove use of BerkeleyDB and DBM::Deep in favor of SQLite, like in Transposome
script: find_graph_edges.pl
 - [ ] update to use latest indexing and edge-finding methods
script: repbase_to_typemap.pl
 - [ ] update to use latest mapping methods; incorporate unknowns
script: pfam_fetch
 - [ ] get GO terms and other metadata with models
    
## phylogenetics
script: eutils_taxonomy_methods.pl
 - [ ] remove use of bioperl and use REST API to eutils
script: orthoMCL...
 - [ ] remove version numbers and clean up orthomcl scripts
script: smith-waterman.pl
 - [ ] initialize diagonals in matrix before incrementing
script: cpbase_fetch.pl
 - [ ] remove options that are not working, such as retrieving the alignments