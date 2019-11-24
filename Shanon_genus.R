#### ASv based method with ASV shannon diversity of non focal genera (from unrarefied dataset)
#### ASV diversity of the focal genus is calculated from rarefied dataset
#### glmm: ASV diversity of the focal genus ~ ASV shannon diversity of non focal genera
#### 31 aout 2018 - NM

if(!require(lme4)){install.packages("lme4")}
require(lme4)

if(!require(readxl)){install.packages("readxl")}
require(readxl)

if(!require(vegan)){install.packages("vegan")}
require (vegan)

if(!require(tidyr)){install.packages("tidyr")}
require (tidyr)

if(!require(ggplot2)){install.packages("ggplot2")}
require (ggplot2)

### list of genera and their ASV counts estimated on the unrarefied dataset
df.ASV=read.table("genus_asv_shanon_final/genus_ASV_list_unrare2.tsv",header=T,sep=",")
df.data=df.ASV[,-1]

test1=df.data
### data.wide is the sample*genus matrix 
data.wide <- spread(data = test1, 
             key = genus_type,
             value = nb_ASV,fill=0)

rownames(data.wide)=data.wide[,1]
data.wide=data.wide[,-1]
test1$shanon=0

#### calcul shannon pour chaque genre = shannon des autres genre - le genre focal
#### df.shannon est la matrice sample*genres non focaux
#### shannon doit avoir une matrice site*esp en entrée

for ( focalG in unique(test1$genus_type)) {
  df.shanon=data.wide[,-which(colnames(data.wide)==focalG)]
  shanon=diversity(df.shanon,index="shannon")
  shanon.df=as.data.frame(shanon)  ### shanon.df = 2000 samples*1shannon_genrefocal
  shanon.df$sample=rownames(shanon.df)
  rownames(shanon.df)=c()
  colnames(shanon.df)=c("shannon","sample")
  j=1
  for (j in 1:dim(shanon.df)[1]) {
    test1[(test1$genus_type == focalG) & (test1$sample == shanon.df[j,"sample"]),"shanon"] = shanon.df[j,"shannon"]
  }
}

### test1 = matrice qui pour chaque genre , elle a autant d'entrées que de samples où il est présent
#### avec le shannon indice des genera non focaux dans le dit sample

### merge test1 avec le dataset rarefié tel qu'on ait dans l'analyse nb_ASV rarefié et shanon non rarefié
df.genus=read.table("genus_ASV_list.tsv",header=T,sep=",")
tab.genus.shanon=test1[,c("sample","shanon","genus_type")]
df.glm=merge(df.genus,test1[,c("sample","shanon","genus_type")],by=c("sample","genus_type"))


################## studies list
study_lab=read_excel("lab_emp.xlsx", col_names  = TRUE)
colnames(study_lab)=c("study","PI")

##### table des empo pour chaque sample
tab_emp=read.csv("sample_emplevel.tsv",header=T)
colnames(tab_emp)=c("sample","empo_0", "empo_1" ,    "empo_2"  ,   "empo_3")

tab_emp=as.data.frame(tab_emp)
tab_emp$empo_1=as.factor(tab_emp$empo_1)
tab_emp$empo_2=as.factor(tab_emp$empo_2)
tab_emp$empo_3=as.factor(tab_emp$empo_3)

tab.data=df.glm
tab.data$genus_type=as.factor(tab.data$genus_type)

tab.data_emp=merge(tab.data,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")

for (i in 1:nrow(tab.data_emp) )
  tab.data_emp[i,'study']=strsplit(as.character(tab.data_emp[i,]$sample),"[.]")[[1]][1]

tab.final=merge(tab.data_emp,study_lab,by="study")
tab.final$PI=as.factor(tab.final$PI)
tab.final$genus_type=as.factor(tab.final$genus_type)
tab.final$empo_3=as.factor(tab.final$empo_3)

## Z standardisation
pvar='shanon'
datsc1 <- tab.final
datsc1[pvar] <- lapply(tab.final[pvar],scale)
head(datsc1)
colnames(tab.final)

### GLMM 
mymodel=glmer(nb_ASV~shanon+(shanon|genus_type/empo_3)+(shanon|empo_3/PI)+(shanon|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(glmer.genus1)
#save(mymodel,file="full_model_genusASV_UnrarShanon.RData")
#load("/Users/naima/Projet_Diversification/ASV_method/genus_asv_shanon_final/full_model_genusASV_UnrarShanon.RData")
summary(mymodel)

##enlever la singularité empo3 (cor=1) et sample (Cor=1) 
mymodel.1=glmer(nb_ASV~shanon+(shanon|genus_type/empo_3)+(shanon|empo_3:PI)+(shanon-1|empo_3)+(shanon-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
save(mymodel.1,file="mymodel_1.RData")
summary(mymodel.1)
save(mymodel.1,file="mymodelpoint1.RData")

##enlever empo3 car son effet=0 (le model le plus significant)
mymodel.1=glmer(nb_ASV~shanon+(shanon|genus_type/empo_3)+(shanon|empo_3:PI)+(shanon-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(mymodel.1)
save(mymodel.1,file="mymodel1sansemp3.RData")

###signif random
mymodel.2=glmer(nb_ASV~shanon+(shanon|genus_type/empo_3)+(shanon|empo_3:PI),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
mymodel.3=glmer(nb_ASV~shanon+(shanon|genus_type/empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
mymodel.4=glmer(nb_ASV~shanon+(shanon|genus_type:empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
mymodel.5=lm(nb_ASV~shanon,datsc1)
anova(mymodel.1,mymodel.2,mymodel.3,mymodel.4,mymodel.5)

### tester fixed effects
mymodel_rand=glmer(nb_ASV~1+(shanon|genus_type/empo_3)+(shanon|empo_3:PI)+(shanon-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
save(mymodel_rand,file="mymodel_rand.RData")
summary(mymodel_rand)
anova(mymodel_rand,mymodel.1)

overdisp_fun(mymodel.1)  

## plot residual contre independent variable
ggplot(data.frame(x1=datsc1$shanon,pearson=residuals(mymodel.1,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw()

##
plot(fitted(mymodel.1), residuals(mymodel.1),xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)

###
qqnorm(scale(resid(mymodel.1)),ylab="Residual quantiles",col="orange")
qqline(scale(resid(mymodel.1)),col="blue")

