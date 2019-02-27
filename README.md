# DBD in the EMP dataset

dbd_analys_input.py : associates for every sample, the taxonomic annotations of all the ASVs it contains. 

glmm_analys_input.py : uses the output of glmm_analys_input.py to estimate diversification rate for every taxonomic ratio (ASV:genus, genus:family, family:order, order:class and class:phylum), as welle as diversity (the number of non-focal lineages: non focal phyla, classes, orders, families and genera numbers).

Construct_fasta_per_sample.py : build the Fasta file for every sample from the EMP summary file. 
