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

###ajouter un booleen pour native/migrant/generalist
for (i in 1: dim(df.animales)[1]) {
  if (df.animales[i,'genus_type'] %in% unique(animales$genus_type) )
      df.animales[i,'status'] = "Native"
  else
    if (df.animales[i,'genus_type'] %in% unique(Saline$genus_type) || df.animales[i,'genus_type'] %in% unique(nonsaline$genus_type))
        df.animales[i,'status'] = "Migrant"
    else
      df.animales[i,'status'] = "Generalist"
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
glmm.animal.1 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3/PI)+(nb_genus|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

# summary(glmer.full.status.1)
 save(glmm.animal.1,file="glmm_animal1.RData")

###enlever  singularité sur sample et empo_3 car cor=1
glmm.animal.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+(nb_genus-1|empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
save(glmm.animal.2,file="glmm_animal2.RData")

###enlever empo3 car 0 (et le plus significant)
glmm.anim.3 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
save(glmm.anim.3,file="glmm_animal3.RData")
summary(glmm.anim.3)


### Random effect significance
#glmm.anim.3 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.4 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.5 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.6 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type:empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))

glmer.full.7 = lm(nb_ASV~nb_genus*status,datsc1)

anova(glmm.anim.3,glmer.full.4,glmer.full.5,glmer.full.6,glmer.full.7)
# Df    AIC    BIC logLik deviance     Chisq Chi Df Pr(>Chisq)    
# glmer.full.7  7 136314 136373 -68150   136300                                
# glmer.full.6  9 100007 100083 -49995    99989 36310.980      2  < 2.2e-16 ***
#   glmer.full.5 12  99341  99442 -49658    99317   672.584      3  < 2.2e-16 ***
#   glmer.full.4 15  98770  98896 -49370    98740   576.606      3  < 2.2e-16 ***
#   glmm.anim.3  16  98718  98853 -49343    98686    53.733      1  2.297e-13 ***

###tester la significance de l'interaction
glmer.full.31 = glmer(nb_ASV~nb_genus+status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(glmm.anim.3,glmer.full.31)
# Df   AIC   BIC logLik deviance Chisq Chi Df Pr(>Chisq)    
# glmer.full.31 14 98771 98889 -49371    98743                             
# glmm.anim.3   16 98718 98853 -49343    98686 56.273      2  6.031e-13 ***

## plot the interaction
require(sjPlot)
pdf("animal3status.pdf")
plot1=plot_model(glmm.anim.3,type="int",terms=c("nb_genus","status"),title="Animal",show.legend = TRUE, axis.title=c("","ASV number/genus"),colors = c("black","red","orange"))
plot1
dev.off()

####
overdisp_fun(glmm.anim.3)
# chisq        ratio          rdf            p 
# 1.756283e+04 5.274597e-01 3.329700e+04 1.000000e+00


## diagnostic plots
plot(fitted(glmm.anim.3),residuals(glmm.anim.3),xlab = "Fitted Values", ylab = "Residuals")
qqnorm(scale(resid(glmm.anim.3)),ylab="Residual quantiles",col="orange")

