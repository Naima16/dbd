# DBD in the EMP dataset

1. Construct_fasta_per_sample.py : build Fasta files from the EMP summary file. 

2. dbd_analys_input.py : associates for every sample, the taxonomic annotations of all the ASVs it contains. 

3. glmm_analys_input.py : uses the output of glmm_analys_input.py to estimate diversification rate for every taxonomic ratio   (ASV:genus, genus:family, family:order, order:class and class:phylum), as welle as diversity (the number of non-focal lineages: non focal phyla, classes, orders, families and genera numbers).


