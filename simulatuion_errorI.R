# Simulations for type 1 error.

#We will permute the spiders data at random [I used method 1: Permute data in the entire matrix, and method 2, Permute data in individual columns], producing each time a data matrix in which H0 is certainly true. 
#We will then run function dbd.perm2() and obtain the associated p-value.
#We will repeat this simulation a large number of times and obtain an estimate of the proportion of the simulations where H0 was rejected at a variety of alpha rejection levels {0.01, 0.05, 0.10, 0.20}.

#“A statistical testing procedure is valid if the probability of a type I error (rejecting H0 when true) is no greater than α, the level of significance, for any α.” (Edgington 1995, p. 37).


simul.type1.dbd = function(nsim=100, n=28, p=12, g=NULL, pg = 6, gen.method=1, prob.mat1=0.8, sd1=1.5,perm.method=2)
  {
    # Generate matrix mat1 (site x species)
    if(gen.method==1) { # Random Poisson data
      mat1 <- matrix(rpois((n*p),prob.mat1),n,p)
    } else {        # Random lognormal data
      mat1 <- matrix(round(exp(rnorm((n*p),0,sd1))),n,p)
    }
    out = matrix(NA, nsim, 3)
    rownames(out) = paste("Simul",1:nsim,sep=".")
    colnames(out) = c("cor.sg", "cor.t", "p.val")
    
    aa <- system.time({
      for(i in 1:nsim) {
        if(perm.method==1) {
          # Permute data in individual columns of mat1 (sites x species)
          mat1.H0 = apply(mat1,2,sample)
        }
        if(perm.method==2) { # Permute p/a data in individual rows of mat1 (sites x species)
          for(ii in 1:n) mat1.H0[ii,] = sample(mat1[ii,])
        }
        
        tmp = dbd.perm2(mat1.H0, g, perm.method=perm.method, nperm=999, clock=FALSE)
        out[i,] = c(tmp$Stat.out[1,1:2], tmp$p.val)
      }
    })
    aa[3] <- sprintf("%2f", aa[3])
    cat("Time for computation =", aa[3], " sec\n")
    #
    data.frame(out)
  }


dbd.perm2 <- function(mat, g, perm.method=1, nperm=999, nn.print=0,clock=TRUE)

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
      # if(perm.method==1) { # Permute p/a data in in the whole matrix mat1 (sites x species)
      #   mat1.perm=.Call("sampleC_all",mat1)
      # }
      ##this is model 1 in the paper
      if(perm.method==1) { # Permute data in individual columns of mat1 (sites x species)
        mat1.perm = .Call("sampleC_col_real",mat1)
      }
      ## this is model 2 in the paper
      if(perm.method==2) { # Permute p/a data in individual rows of mat1 (sites x species)
        mat1.perm=.Call("sampleC_row_real",mat1)
      }
      
      mat1.perm=.Call("presence_absence",mat1.perm)
      # Create matrix "mat2" of sites x genera
      mat2=.Call("produit",mat1.perm,A) ##produit matriciel
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
  list(Richness=Richness.emp, Stat.out=Stat.out, p.val=p.val,Richness.perm=Richness)
}


genera = c("Alop", "Alop", "Alop", "Arct", "Arct", "Aulo", "Pard", "Pard", "Pard", "Pard", "Troc", "Zora")

res1 = simul.type1.dbd(nsim=100, g=genera, gen.method=2,perm.method=2)

head(res1)
res1$p.val
length(which(res1$p.val <= 0.05))
length(which(res1$p.val <= 0.10))
length(which(res1$p.val <= 0.20))
length(which(res1$p.val <= 0.50))
mean(res1$p.val)



