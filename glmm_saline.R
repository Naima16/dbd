#### glmm ASV:genus pour saline group
#### resident-migrant analysis
#### 7 avril 2019


library(lme4)
library(ggpmisc)
library(calibrate)
library("readxl")
library(car)
library(ggplot2)
library(labdsv)

overdisp_fun <- function(model) {
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

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

######################
##########

#### garder seulement les samples du groupe saline
df.final=tab.final[which(tab.final$empo_3  %in% c("Surface (saline)","Water (saline)","Sediment (saline)","Hypersaline (saline)","Plant surface" )),]


################## les genera residents de salines  dans les 17 empo_3
Saline=read.table("~/genus_ASV_SalineGrp.tsv",header=T,sep=",")

df.nativ.emp=merge(Saline[,-1],tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")

df.native=df.nativ.emp[which(df.nativ.emp$empo_3 %in% c("Surface (saline)","Water (saline)","Sediment (saline)","Hypersaline (saline)","Plant surface" )),]

as.data.frame(table(df.native$empo_3))
df.native.final=unique(df.native[,c("sample","nb_genus")])


#### prendre le #genus des resident et #ASV resident+migrants

df.salines=merge(df.final[,-4],df.native.final,by="sample")

###ajouter un booleen pour resident/migrant/generalist
for (i in 1: dim(df.salines)[1]) {
  if (df.salines[i,'genus_type'] %in% unique(Saline$genus_type) )
    df.salines[i,'status'] = "Native"
  else
    if (df.salines[i,'genus_type'] %in% unique(nonsaline$genus_type) || df.salines[i,'genus_type'] %in% unique(animal$genus_type) )
        df.salines[i,'status'] = "Migrant"
    else
      df.salines[i,'status'] = "Generalist"
}

df.salines$status=as.factor(df.salines$status)

################
############ standardise predictor
pvar='nb_genus'
datsc1 <- df.salines
datsc1[pvar] <- lapply(df.salines[pvar],scale)

datsc_saline=datsc1
save(datsc_saline,file="datsc_saline.RData")

#### full gLMM
glmm.sal.1 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3/PI)+(nb_genus|sample),datsc1,family=poisson(link=log),
                           control=glmerControl(optimizer="bobyqa"))
summary(glmm.sal.1)
save(glmm.sal.1,file="glmm_saline_1.RData")


## enlever la singularité sur sample et empo_3 (le plus signiifcant)
glmm.sal.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+ (nb_genus-1|empo_3)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                control=glmerControl(optimizer="bobyqa"))

summary(glmm.sal.2)
save(glmm.sal.2,file="glmm_saline_2.RData")

###supprimer empo_3 car var 0
glm.status3.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                control=glmerControl(optimizer="bobyqa"))
save(glm.status3.2,file="glm.saline_3.RData")
summary(glm.status3.2)

##tester la significance des random d'après ben bolker : https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
#glm.status3.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
  #                     control=glmerControl(optimizer="bobyqa"))

glm1=update(glm.status3.2,.~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI))
glm2=update(glm.status3.2,.~nb_genus*status+(nb_genus|genus_type/empo_3))
glm3 <- lm(nb_ASV~nb_genus*status,datsc1)
anova(glm.status3.2 ,glm1,glm2,glm3)


# tester l'interaction:
glmer.inter = glmer(nb_ASV~nb_genus+status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                    control=glmerControl(optimizer="bobyqa"))
anova(glm.status3.2,glmer.inter)
# Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# glmer.inter   14 57123 57233 -28547    57095                             
# glm.status3.2 16 57100 57226 -28534    57068 26.994      2  1.375e-06 ***

overdisp_fun(glm.status3.2)
# chisq        ratio          rdf            p 
# 1.029538e+04 5.282932e-01 1.948800e+04 1.000000e+00 

plot(fitted(glm.status3.2),residuals(glm.status3.2),xlab = "Fitted Values", ylab = "Residuals")
qqnorm(scale(resid(glm.status3.2)),ylab="Residual quantiles",col="orange")

##interacton plot
require(sjPlot)
pdf("saline.pdf")
plot1=plot_model(glm.status3.2,type="int",terms=c("nb_genus","status"),title="Saline",show.legend = TRUE, axis.title=c("","ASV number/genus"),colors = c("black","red"))
plot1
dev.off()
