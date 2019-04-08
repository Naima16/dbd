#### glmm resident versus migrant analysis
#### animal cluster 
#### NM 7 avril 2019


library(lme4)
library(ggpmisc)
library(calibrate)
library("readxl")
library(car)
library(ggplot2)
library(labdsv)


### table des PI
study_lab=read_excel("~/lab_emp.xlsx", col_names  = TRUE)
colnames(study_lab)=c("study","PI")


##### table des empo pour chaque sample
tab_emp=read.csv("~/sample_emplevel.tsv",header=T)
tab_emp=tab_emp[,-1]
colnames(tab_emp)=c("sample","empo_0", "empo_1" ,    "empo_2"  ,   "empo_3")

tab_emp=as.data.frame(tab_emp)
tab_emp$empo_1=as.factor(tab_emp$empo_1)
tab_emp$empo_2=as.factor(tab_emp$empo_2)
tab_emp$empo_3=as.factor(tab_emp$empo_3)

##table des ASV-genus
tab.data=read.table("~/genus_ASV_list.tsv",header=T,sep=",")
tab.data=as.data.frame(tab.data)

tab.data_emp=merge(tab.data[,-1],tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")

###ajouter le lab
for (i in 1:nrow(tab.data_emp) )
  tab.data_emp[i,'study']=strsplit(as.character(tab.data_emp[i,]$sample),"[.]")[[1]][1]

tab.final=merge(tab.data_emp,study_lab,by="study")

tab.final$PI=as.factor(tab.final$PI)
tab.final$genus_type=as.factor(tab.final$genus_type)
tab.final$empo_3=as.factor(tab.final$empo_3)


##########

#### garder seulement les samples des  biomes du groupe animal resultats fuzzy k-mean clustering
df.final=tab.final[which(tab.final$empo_3  %in% c("Animal distal gut","Animal proximal gut","Animal secretion","Animal corpus","Animal surface","Aerosol (non-saline)","Plant corpus" )),]


################## les genera natives de animales dans les 17 empo_3
animales=read.table("~/genus_ASV_AnimalGrp.tsv",header=T,sep=",")
df.nativ.emp=merge(animales[,-1],tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")
df.native=df.nativ.emp[which(df.nativ.emp$empo_3 %in% c("Animal distal gut","Animal proximal gut","Animal secretion","Animal corpus","Animal surface","Aerosol (non-saline)","Plant corpus" )),]
as.data.frame(table(df.native$empo_3))

##diversité en native genera pour tous les samples du cluster animal
df.native.final=unique(df.native[,c("sample","nb_genus")])

#### prendre la diversité en genus des natives et le nombre d'ASV native+migrate
#### colonne 4 = nb_genus
df.animales=merge(df.final[,-4],df.native.final,by="sample")

###ajouter un booleen pour native versus migrants
for (i in 1: dim(df.animales)[1]) {
  if (df.animales[i,'genus_type'] %in% unique(animales$genus_type) )
      df.animales[i,'status'] = "Native"
  else
      df.animales[i,'status'] = "Migrant"
}

df.animales$status=as.factor(df.animales$status)

################
############ standardisation (mean=0 and sd=1)
pvar='nb_genus'
datsc1 <- df.animales
datsc1[pvar] <- lapply(df.animales[pvar],scale)

datsc_animal=datsc1
save(datsc_animal,file="datsc_animal.RData")


########################### Full glmm
glmer.full.status.1 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3/PI)+(nb_genus|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

summary(glmer.full.status.1)
save(glmer.full.status.1,file="glmm_animal1.RData")

###enlever  singularité sur sample et empo_3 car cor=1
glmer.full.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+(nb_genus-1|empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
save(glmer.full.2,file="glmm_animal2.RData")

###enlever empo3 car 0 (et le plus significant)
glmer.full.3 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
save(glmer.full.3,file="glmm_animal3.RData")
summary(glmer.full.3)


### random effect significance
glmer.full.4 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.5 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.6 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type:empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.7 = lm(nb_ASV~nb_genus*status,datsc1)

anova(glmer.full.3,glmer.full.4,glmer.full.5,glmer.full.6,glmer.full.7)
#   glmer.full.7  5 136317 136359 -68153   136307                                
#   glmer.full.6  7 100005 100064 -49996    99991 36315.902      2  < 2.2e-16 ***
#   glmer.full.5 10  99338  99423 -49659    99318   672.572      3  < 2.2e-16 ***
#   glmer.full.4 13  98770  98879 -49372    98744   574.914      3  < 2.2e-16 ***
#   glmer.full.3 14  98717  98835 -49345    98689    54.386      1  1.647e-13 ***

###tester la significance de l'interaction
glmer.full.31 = glmer(nb_ASV~nb_genus+status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(glmer.full.31,glmer.full.3)
# Df   AIC   BIC logLik deviance Chisq Chi Df Pr(>Chisq)    
# glmer.full.31 13 98771 98881 -49373    98745                            
# glmer.full.3  14 98717 98835 -49345    98689 56.01      1   7.21e-14 ***

require(sjPlot)
pdf("animal.pdf")
plot1=plot_model(glmer.full.3,type="int",terms=c("nb_genus","status"),title="Animal",show.legend = TRUE, axis.title=c("","ASV number/genus"),colors = c("black","red"))
plot1
dev.off()

# ###
overdisp_fun(glmer.full.3)
# chisq        ratio          rdf            p 
# 1.756283e+04 5.274597e-01 3.329700e+04 1.000000e+00

## diagnostic plots
plot(fitted(glmer.full.3),residuals(glmer.full.3),xlab = "Fitted Values", ylab = "Residuals")
qqnorm(scale(resid(glmer.full.3)),ylab="Residual quantiles",col="orange")

