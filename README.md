# DBD in the EMP dataset

1. Construct_fasta_per_sample.py : builds Fasta files from the EMP summary file. 

2. dbd_analys_input.py : for every sample, it associates the taxonomic annotations of all the ASVs it contains. 

3. glmm_analys_input.py : uses the output of dbd_analys_input.py to estimate diversification rate for every taxonomic ratio   (ASV:Genus, Genus:Family, Family:Order, Order:Class and Class:Phylum), as well as diversity (the number of non-focal lineages: non focal phyla, classes, orders, families and genera numbers).

4. Figure 3 : does the indicator species analysis and the PCA, then plot 1) figure 3-A for visualizing the clustering and the indicator genera for every cluster; 2) figure 3-B using GLMM output (a GLMM for every environment cluster : animal-associated, saline and non-saline).

4. percentidentity.py : clusters ASVs by different percent identity using USEARCH.

5. Figure4.r : plot the figure 4 using GLMMs output (1: glmm with the interaction between diversity and biome as fixed effect to predict diversification (supplementary, file1, section 5), and 2: the GLMM of the genome size analysis (supplementary, file1, section 5)

5. GLMM_CLASSperPHYLUM.R : Generalized linear mixed model for Class:Phylum ratio for diversity begets diversitification  analysis (Table 1) as well as DBD slope variation across different biomes analysis (Table 4-A).


