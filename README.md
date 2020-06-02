# DBD in the EMP dataset

1. summarize_otu_distributions.py : generates  otu_summary.tsv from biome table (this script is from Thompson et al. 2017).

2. dbd_analys_input.py : for every sample, it associates the taxonomic annotations of all the ASVs it contains.  

3. glmm_analys_input.py : uses the output of dbd_analys_input.py to estimate diversification rate for every taxonomic ratio   (ASV:Genus, Genus:Family, Family:Order, Order:Class and Class:Phylum), as well as diversity (the number of non-focal lineages: non focal phyla, classes, orders, families and genera numbers).

4. GLMM_CLASSperPHYLUM.R : Generalized linear mixed model for diversity begets diversitification  analysis (Table 1) as well as DBD slope variation across different biomes analysis (Supplementary Data file 1 Section 4.3), for Class:Phylum ratio. For the other ratio, the table in line 47 have to be replaced by the appropriate diversity_diversification table (glmm_analys_input.py output).

5. Construct_fasta_per_sample.py : builds Fasta files from the EMP summary file. 

6. Figure 3 : Indicator species analysis and PCA, and plot 1) figure 3-A for visualizing the clustering and the indicator genera for every cluster; 2) figure 3-B using GLMM output (a GLMM for every environment cluster : animal-associated, saline and non-saline).

7. percentidentity.py : clusters ASVs by different percent identity using USEARCH.

8. Figure4.r : plot the figure 4 using GLMMs output (1: glmm with the interaction between diversity and biome as fixed effect to predict diversification (supplementary, file1, section 5), and 2: the GLMM of the genome size analysis (supplementary, file1, section 5)

9. GLMM_animal.R : Generalized linear mixed model for ASV:Genus in animal cluster (resident-migrant analysis).

10. GLMM_saline.R : Generalized linear mixed model for ASV:Genus in saline cluster (resident-migrant analysis).

11. GLMM_NonSaline.R : Generalized linear mixed model for ASV:Genus in non-saline cluster (resident-migrant analysis).

12. identity_permtest.R : permutation test for diversification~diversity relationship in Nc percent identity analysis.

13. FigureS10.R : plot Figure S10, linear, quadratic and cubic models for diversification~diversity relationship (Nc percent identity analysis)

14. glmm_abiotic_phylum_class.r : Abiotic factors GLMM for Class:Phylum ratio. This script may be used for any ratio, if the table with phyla count (diversity) and class:phyla ratio (diversification) for every phylum in all samples (line 33) is replaced (by diversity-diversification table for any other ratio (glmm_analys_input.py output)).

15. glmm_GenomeSize.r : Genome size analysis.
16. permute_ASVs_synthetic.pl : Rarefaction simulation.


