library(ggplot2)

###spiders. test example
spiders = read.table("/Users/naima/DBD_permutation_EMP/16sept/script_finaux_18sept/biom_to_csv/25sept/Spiders_28x12_spe.txt")
genera = c("Alop", "Alop", "Alop", "Arct", "Arct", "Aulo", "Pard", "Pard", "Pard", "Pard", "Troc", "Zora")
mat = as.matrix(spiders)
genus=genera

### emp
EMP_data=read.csv("/Users/naima/Projet_Diversification/Dernier_nov/permutation_9oct_2019/DBD_permutation_EMP/16sept/script_finaux_18sept/biom_to_csv/EMP_table.from_biom_w_taxonomy.txt",header=T,sep="\t")
dim(EMP_data)

emp1=EMP_data
emp1=emp1[,-1]

genus=array()
for (i in 1:dim(emp1)[1])
  genus[i]=unlist(strsplit(unlist(strsplit(toString(emp1[i,"taxonomy"]),";"))[6],"__"))[2]

nb_genus=length(genus[!is.na(genus)]) 

emp2=cbind(emp1,genus)
#### drop unnamed genera 
emp2=emp2[ !is.na(emp2$genus),]
colnames(emp2)=c(colnames(emp1),"genus")

ind_taxonomy=which(colnames(emp2)=="taxonomy")
ind_genus=which(colnames(emp2)=="genus")

mat = t(emp2[,-c(ind_taxonomy,ind_genus)])
genus=emp2$genus

### C code import
dyn.load("/Users/naima/DBD_permutation_EMP/1oct_resultfinaux_alltax/25sept_genus_Asv/sample2.so")


nperm=99

out=dbd.perm2(mat,g=genus,perm.method=2,nperm=nperm,nn.print=0,clock=T)
save(out,file="gen_asv_perm_method2.RData")

out=dbd.perm2(mat,g=genus,perm.method=3,nperm=nperm,nn.print=0,clock=T)
save(out,file="genus_asv_perm_model3.RData")


###
dbd.perm2 <- function(mat, g, perm.method=1, nperm=999, nn.print=0,clock=TRUE)
  # Permutation test for the "Diversity begets diversity" (dbd) paper.
  # Vers. 2: general permutation of values by rows, by columns, or in the whole matrix "mat". Statistic to be tested: rich.genera ~ species.richness/N.genera 
  # 
  # Arguments –
  # mat = data matrix (data.frame or matrix), sites in rows x species (or OTUs) in columns.
  # g = vector of {genera, families, orders, classes, phyla} names in columns of mat
  # perm.method = {1,2,3} : Different permutation methods. 
  # nperm = number of permutations.
  # nn.print : print 'nn' lines of the output matrix. Default: nn=0.
  # clock : If \code{clock=TRUE}, the computation time is printed. Useful when nperm large.
# 
# Details –
# Matrix "mat2" (sites x genera) is produced first for the unpermuted data (mat2 = mat %*% A.perm), then after each permutation using one of the three methods below. 
# 
# Statistic to be tested: richness.genera ~ (how many species per genus in each site).
#
# A permutation test is a valid way of testing the relationship between x and y/x. The permutations must be done on y, or on the separate elements y and x, before the ratio variable y/x is recomputed after each permutation. Reference: Jackson & Somers (1991).
#
# Three permutation methods are implemented –
# Before computing mat2,
# Method 1. Permute p/a data in in the whole matrix mat1. (dropped from the paper)
# Method 2. Permute data in individual columns of mat1. (Model 1 in the paper)
# Method 3. Permute data in individual rows of mat1. (Model 2 in the paper)
#    ### Method 2 has low power. ###

# Value –
# The output list comprises the following elements:
# • Richness : data frame with 2 columns, rich.gen and Nsp.by.Ngen, for data.
# • Stat.out contains 2 columns: the correlation and the associated t-statistic.
# • p.val is the p-value computed from the test of the correlation statistic. 
# • Richness.perm : data frame with 2 columns, rich.gen and Nsp.by.Ngen, permuted data.
# 
# References –
# Jackson, D. A. & K. M. Somers. 1991. The spectre of ‘spurious’ correlations. 
# Oecologia 86: 147-151.
#
# P. Legendre, Naima Madi september 2019
{
  epsilon <- sqrt(.Machine$double.eps)
  aa <- system.time({
    mat = as.matrix(mat)
    n <- nrow(mat)
    p <- ncol(mat)
    # Create matrix mat1 (sites x species), p/a data
    mat1=.Call("presence_absence",mat)
    rich.sp=.Call("richesse",mat1)
    # Create matrix A of dummy variables (species x genera)
    g = as.factor(g)
    ng = length(g)
    A = model.matrix(~g-1)
    A = A[1:ng,]
    colnames(A) = levels(g)
    # Create matrix mat2 (sites x genera), p/a data
    mat2=.Call("produit",mat,A)
    mat2=.Call("presence_absence",mat2)
    rich.gen=.Call("richesse",mat2)
     
    Nsp.by.Ngen = rich.sp/rich.gen
    Richness = as.data.frame(cbind(rich.gen,rich.sp, Nsp.by.Ngen))
    Richness=Richness[Richness$rich.gen !=0,]
    Richness.emp=Richness
    Stat.out = matrix(NA,nperm+1,2)
    colnames(Stat.out) = c("cor.sg", "cor.t")
    rownames(Stat.out) = paste("perm", 0:nperm, sep=".")
    
    # Compute statistics for unpermuted data and write them to Stat.out
    Stat.out[1,1] <- cor.sg <- cor(Richness$rich.gen, Richness$Nsp.by.Ngen) 
    Stat.out[1,2] <- cor.t  <- cor.sg * sqrt(n-2)/sqrt(1-cor.sg^2)
    if(cor.sg > 1.0-sqrt(epsilon)) cat("cor = 1.0; a test is useless.\n")
    #
    # Permutations begin
    if((perm.method<1) | (perm.method>3)) stop ("Method not implemented")
    for(i in 1:(nperm)) {
      
      mat1.perm = matrix(0,n,p)
      ### this permutation was dropped from the paper
      if(perm.method==1) { # Permute p/a data in in the whole matrix mat1 (sites x species)
        mat1.perm=.Call("sampleC_all",mat1)
      }
      ##this is model 1 in the paper
      if(perm.method==2) { # Permute data in individual columns of mat1 (sites x species)
        mat1.perm = .Call("sampleC_col_real",mat1)
      }
      ## this is model 2 in the paper
      if(perm.method==3) { # Permute p/a data in individual rows of mat1 (sites x species)
        mat1.perm=.Call("sampleC_row_real",mat1)
      }
  
      # Create matrix "mat2" of sites x genera
      mat2=.Call("produit",mat1.perm,A) 
      mat2=.Call("presence_absence",mat2)
      rich.gen=.Call("richesse",mat2)
      length(rich.gen[rich.gen==0])
      rich.sp=.Call("richesse",mat1.perm)
      
      Nsp.by.Ngen.perm = rich.sp/rich.gen
      
      Richness = as.data.frame(cbind(rich.gen, Nsp.by.Ngen.perm))
      
      Richness=Richness[Richness$rich.gen !=0,]
      # Compute statistics
      Stat.out[i+1,1] <- cor.sg.perm <- cor(Richness$rich.gen, Richness$Nsp.by.Ngen.perm) 
      Stat.out[i+1,2] <- cor.t.perm  <- cor.sg.perm * sqrt(n-2)/sqrt(1-cor.sg.perm^2)
    }
    #
    # Print a few lines (nn lines) of the output matrix
    if(nn.print>0) print(Stat.out[1:nn.print,])
  })
  aa[3] <- sprintf("%2f", aa[3])
  if (clock) cat("Time for computation =", aa[3], " sec\n")
  
  # Compute p-value
  p.val = length(which(Stat.out[,2] >= Stat.out[1,2]))/(nperm+1)
  #
  ###on devrait sortir le Richness des vraies data celui la est permuté
  list(Richness=Richness.emp,Stat.out=Stat.out, p.val=p.val,Richness.perm=Richness)
}
