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

###ajouter un booleen pour resident versus migrants
for (i in 1: dim(df.salines)[1]) {
  if (df.salines[i,'genus_type'] %in% unique(Saline$genus_type) )
    df.salines[i,'status'] = "Native"
  else
    df.salines[i,'status'] = "Migrant"
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
glmer.full.status.Sal.1 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3/PI)+(nb_genus|sample),datsc1,family=poisson(link=log),
                           control=glmerControl(optimizer="bobyqa"))
summary(glmer.full.status.Sal.1)
save(glmer.full.status.Sal.1,file="glmm_saline1.RData")


## enlever la singularité sur sample et empo_3 (le plus signiifcant)
glmer.full.status.Sal.21 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+ (nb_genus-1|empo_3)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                control=glmerControl(optimizer="bobyqa"))

summary(glmer.full.status.Sal.21)
save(glmer.full.status.Sal.21,file="glmm_saline2.RData")

###supprimer empo_3 car var 0
glmer.full.status.Sal.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                control=glmerControl(optimizer="bobyqa"))
save(glmer.full.status.Sal.2,file="glmer.full.status.Sal.signif.RData")


##tester la significance des random d'après ben bolker : https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
glmer.full.status.Sal.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                control=glmerControl(optimizer="bobyqa"))

glm1=update(glmer.full.status.Sal.2,.~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI))
glm2=update(glmer.full.status.Sal.2,.~nb_genus*status+(nb_genus|genus_type/empo_3))
glm3 <- lm(nb_ASV~nb_genus*status,datsc1)
anova(glmer.full.status.Sal.2 ,glm1,glm2,glm3)

# Df   AIC   BIC logLik deviance     Chisq Chi Df Pr(>Chisq)    
# glm3                     5 77116 77155 -38553    77106                                
# glm2                    10 57313 57392 -28647    57293 19812.712      5  < 2.2e-16 ***
# glm1                    13 57136 57239 -28555    57110   183.021      3  < 2.2e-16 ***
# glmer.full.status.Sal.2 14 57106 57217 -28539    57078    31.823      1  1.689e-08 ***

# interaction:
glmer.inter = glmer(nb_ASV~nb_genus+status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                control=glmerControl(optimizer="bobyqa"))
anova(glmer.full.status.Sal.2,glmer.inter)

overdisp_fun(glmer.full.status.Sal.2)
# chisq        ratio          rdf            p 
# 1.029538e+04 5.282932e-01 1.948800e+04 1.000000e+00 

plot(fitted(glmer.full.status.Sal.2),residuals(glmer.full.status.Sal.2),xlab = "Fitted Values", ylab = "Residuals")
qqnorm(scale(resid(glmer.full.status.Sal.2)),ylab="Residual quantiles",col="orange")

##interacton plot
require(sjPlot)
pdf("saline.pdf")
plot1=plot_model(glmer.full.status.Sal.2,type="int",terms=c("nb_genus","status"),title="Saline",show.legend = TRUE, axis.title=c("","ASV number/genus"),colors = c("black","red"))
plot1
dev.off()
