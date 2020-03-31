##############################################################################################################

#                                       Plots for Naima's paper

#  Baseline dataset with POSITIVE slope
#  run on cluster

##############################################################################################################

#### cluster 

#### load Naima's data
#----------------------
setwd("C:/Users/carme/Dropbox/Work_Files/shapiro_lab/Naima_project/cluster") #laptop
#setwd("/home/cmurall/naimaproject") #cluster

#EMP_data <-read.delim("EMP_table.from_biom_w_taxonomy.txt", header= TRUE) #takes ~3.6 mins
load(file="emp2.RData") 
#or run Naima's code (see "randomGPO_data.R" file for this)


#### libraries
#-------------
library(tictoc)
library(gridExtra)
library(ggplot2)
library(lattice)
library(grid)
library(plyr)
library(dplyr)
library(magrittr)
library(gamlss)


#### used across code
#--------------------
nRows <-22014 # no. of ASVs
nCols <-2000   # no. of samples
num<-nRows*nCols #total no. of elements

p <- c(1, 5504, 11007, 16510, 22014) #c(1, 5504, 11007, 16510, 22014) #(0,25%,50%,75%, 100%) elements are tested  
no<- length(p) #number of plots in a row
pl <-list() #for plots
perlist<-c(0, 25, 50, 75, 100) #for labelling top row of images and images

#### 1. add DBD ----------------------------------------------------------------------------------------------

#create Poisson data, change lambda per site (for groups of columns):
set.seed(3579)
simMat<-matrix(NA,nRows,nCols) #mat for filling
max<-0.01 #max lambda
sites<-read.csv(file="site_structure.csv", header = TRUE) #load Naima's site structure
colnames(sites) <- c("sites", "num_samples") #remove weird character
noSites<-as.vector(sites[,2])
cum<-cumsum(noSites)

for(i in 1:cum[1]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[1]+1):cum[2]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[2]+1):cum[3]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[3]+1):cum[4]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[4]+1):cum[5]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[5]+1):cum[6]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[6]+1):cum[7]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[7]+1):cum[8]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[8]+1):cum[9]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[9]+1):cum[10]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[10]+1):cum[11]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[11]+1):cum[12]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[12]+1):cum[13]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[13]+1):cum[14]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[14]+1):cum[15]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[15]+1):cum[16]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[16]+1):cum[17]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}

write.csv(simMat, file="simMat2_DBD.csv")

# loop to add DBD

rm(i)
for (i in 1:no){
  #ensure still using same base dataset:
  simMat <-read.csv(file="simMat2_DBD.csv", header = TRUE) #or fromGPO
  simMat<-simMat[,-1] #drop first column because it's not data
  
  #make simData a P/A matrix
  simData<-data.frame(simMat) #make into dataframe
  fun<-function(x) replace(x, x>0, 1) #function, replace nonzero elements into 1
  paData <-simData %>% mutate_all(fun) #make data presence/absence
  
  #for filling mulitple elements:
  for (k in 1:nCols){
    for(j in 1:p[i]){ #check p elements per column
      ind<-which(paData[,k]==0, arr.ind = TRUE) #indices of zero elements in column
      if (any(ind!=0)) rind<-sample(ind,1) else next #randomly select a zero element to replace
      el<-sample(paData[,k],1) #choose an element in a column at random
      if (el>0) paData[,k][rind]<-1 #if el is nonzero, fill rind with 1
    }}
  #paData is now transformed
  #save dataset
  percent<-perlist[[i]] #percent label
  name<-paste0("paData2_DBD", percent,".csv") #create label
  write.csv(paData, file=name) #save dataset
  
  #get number of ASVs per sample (by column)
  vecSum <-c()
  for(l in 1:nCols){
    x<-sum(paData[,l])
    vecSum<-append(vecSum,x)
  }
  
  # new dataframe for results
  df <-data.frame()
  df<-cbind(vecSum)
  colnames(df)<-c("no.species")
  
  #get number of genera per sample
  x1<-c()
  for(m in 1:nCols){
    ind <-which(paData[,m]>0, arr.ind = TRUE) #indices for nonzero elements
    vecGen <-emp2$genus[ind] #genera for those indices
    count <-nrow(plyr::count(vecGen))#freq of genera in vecGen, count no. rows which gives no. of genera in vecGen
    x1<-append(x1,count)#store
  }
  x1 <-as.data.frame(x1)
  
  #fill results df
  df <-as.data.frame(df)
  df <- bind_cols(df, x1) #bind
  colnames(df)<-c("no.species","no.gen")
  ratio<- df$no.species/df$no.gen
  df$ratio <-ratio
  head(df)
  
  #plot ASV:#genera vs. # of genera
  pl[[i]]<-ggplot(df, aes(x = no.gen, y=ratio))+geom_point()+
    #xlim(0,12)+#ylim(-5,50)+
    geom_smooth(method = "lm")  + 
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(), plot.title = element_text(hjust = 0.5, face = "plain")) + 
    ggtitle(as.character(percent))
}




#### 2. add EC ----------------------------------------------------------------------------------------------

set.seed(3579)
simMat<-matrix(NA,nRows,nCols) #mat for filling
max<-0.01 #max lambda
#sites<-read.csv(file="site_structure.csv", header = TRUE) #load Naima's site structure
#colnames(sites) <- c("sites", "num_samples")
#noSites<-as.vector(sites[,2])
#cum<-cumsum(noSites)

rm(i) #clear i
for(i in 1:cum[1]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[1]+1):cum[2]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[2]+1):cum[3]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[3]+1):cum[4]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[4]+1):cum[5]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[5]+1):cum[6]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[6]+1):cum[7]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[7]+1):cum[8]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[8]+1):cum[9]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[9]+1):cum[10]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[10]+1):cum[11]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[11]+1):cum[12]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[12]+1):cum[13]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[13]+1):cum[14]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[14]+1):cum[15]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[15]+1):cum[16]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}
for(i in (cum[16]+1):cum[17]){
  lambda<-runif(1, 0, max) #choose a lambda per site (a group of columns)
  simMat[,i]<-rpois(nRows, lambda) #fill elements, each column has its unique distribution 
}

write.csv(simMat, file="simMat2_EC.csv")


# loop to add EC

rm(i)
for (i in 1:no){
  #ensure still using same base dataset:
  simMat <-read.csv(file="simMat2_EC.csv", header = TRUE) #or fromGPO
  simMat<-simMat[,-1] #drop first column because it's not data
  
  #make simData a P/A matrix
  simData<-data.frame(simMat) #make into dataframe
  fun<-function(x) replace(x, x>0, 1) #function, replace nonzero elements into 1
  paData <-simData %>% mutate_all(fun) #make data presence/absence
  
  #for removing mulitple elements:
  for (k in 1:nCols){
    for(j in 1:p[i]){ #check p elements per column
      ind<-which(paData[,k]!=0, arr.ind = TRUE) #indices of nonzero elements in column
      if (any(ind!=0)) rind<-sample(ind,1) else next #randomly select a nonzero element to replace #if (any(ind!=0)) makes sure that ind isn't empty
      el<-sample(paData[,k],1) #choose an element in the column at random
      if (el>0) paData[,k][rind]<-0 else next #if el is nonzero, replace element at position rind with 0
    }}
  #paData is now transformed
  #save dataset
  percent<-perlist[[i]] #percent label
  name<-paste0("paData2_EC", percent,".csv") #create label
  write.csv(paData, file=name) #save dataset
  
  
  #get number of ASVs per sample (by column)
  vecSum <-c()
  for(l in 1:nCols){
    x<-sum(paData[,l])
    vecSum<-append(vecSum,x)
  }
  
  # new dataframe for results
  df <-data.frame()
  df<-cbind(vecSum)
  colnames(df)<-c("no.species")
  
  #get number of genera per sample
  x1<-c()
  for(m in 1:nCols){
    ind <-which(paData[,m]>0, arr.ind = TRUE) #indices for nonzero elements
    vecGen <-emp2$genus[ind] #genera for those indices
    count <-nrow(plyr::count(vecGen))#freq of genera in vecGen, count no. rows which gives no. of genera in vecGen
    x1<-append(x1,count)#store
  }
  x1 <-as.data.frame(x1)
  
  #fill results df
  df <-as.data.frame(df)
  df <- bind_cols(df, x1) #bind
  colnames(df)<-c("no.species","no.gen")
  ratio<- df$no.species/df$no.gen
  df$ratio <-ratio
  head(df)
  
  #plot ASV:#genera vs. # of genera
  pl[[(i+no)]]<-ggplot(df, aes(x = no.gen, y=ratio))+geom_point()+
    #xlim(0,12)+#ylim(-5,50)+
    geom_smooth(method = "lm") + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) #+ ggtitle("rpois,  75% ")
}


#### 3. final plot ------------------------------------------------------------------------------------------

# plot all runs in one plot
fullpl<-grid.arrange(grobs=pl, nrow=2, top="Poisson dataset, by site lambda (0, 0.01)", 
                     left = textGrob("ASV:Genus", rot = 90, vjust = 0.5), 
                     bottom = textGrob("number of genera"))#,
# missing label of rows, DBD top, EC bottom

#save plot
ggsave(fullpl,filename="multiplot2.pdf")#, width = 14.0, height = 5.0, units = "cm")
ggsave(fullpl,filename="multiplot2png.png")
#save data of plot (so it can be reformated)
save(pl, file="multiplot2_data.RData") 
#load using:  load("multiplot2_data.RData")
#then, change like this:
# pl + ggtitle("better title")

# #reconfiguring the multi-plot figure:
# pl[[1]]<-pl[[1]]+ylim(1,2.5)+xlim(0,320)
# pl[[2]]<-pl[[2]]+ylim(1,2.5)+xlim(0,320)
# pl[[3]]<-pl[[3]]+ylim(1,2.5)+xlim(0,320)
# pl[[4]]<-pl[[4]]+ylim(1,2.5)+xlim(0,320)
# pl[[5]]<-pl[[5]]+ylim(1,2.5)+xlim(0,320)
# pl[[6]]<-pl[[6]]+ylim(1,2.)+xlim(0,320)
# pl[[7]]<-pl[[7]]+ylim(1,2.)+xlim(0,320)
# pl[[8]]<-pl[[8]]+ylim(1,2.)+xlim(0,320)
# pl[[9]]<-pl[[9]]+ylim(1,2.)+xlim(0,320)
# pl[[10]]<-pl[[10]]+ylim(1,2.)+xlim(0,320)
