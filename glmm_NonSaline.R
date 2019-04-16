####glmer ASV:genus pour non saline cluster
##  glmer avec genome size 7 avril


library(lme4)
library("readxl")
library(ggplot2)

### table des PI
study_lab=read_excel("~/lab_emp.xlsx", col_names  = TRUE)
colnames(study_lab)=c("study","PI")

##### table des empo pour chaque sample
tab_emp=read.csv("/Users/naima/Projet_Diversification/Diversification_rate_11mars/sample_emplevel.tsv",header=T)
as.data.frame(tab.emp3)
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


#### garder seulement les samples du groupe Non saline 
df.final=tab.final[which(tab.final$empo_3 %in% c("Soil (non-saline)","Surface (non-saline)", "Sediment (non-saline)" ,"Plant rhizosphere" ,"Water (non-saline)")),]

################## les genera residents de NonSaline 
NonSaline=read.table("~/genus_ASV_NonSalineGrp.tsv",header=T,sep=",")

df.nativ.emp=merge(NonSaline[,-1],tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")

df.native=df.nativ.emp[which(df.nativ.emp$empo_3 %in% c("Soil (non-saline)","Surface (non-saline)", "Sediment (non-saline)" ,"Plant rhizosphere" ,"Water (non-saline)")),]
as.data.frame(table(df.native$empo_3))

####garder le nombre de genera resident total pour chaque sample
df.native.final=unique(df.native[,c("sample","nb_genus")])

#### prendre le #genus des residents et #ASV resident+migrate
df.NonSaline=merge(df.final[,-4],df.native.final,by="sample")

###ajouter un booleen pour resident versus migrants
for (i in 1: dim(df.NonSaline)[1]) {
  if (df.NonSaline[i,'genus_type'] %in% unique(NonSaline$genus_type) )
    df.NonSaline[i,'status'] = "Native"
  else
      if (df.NonSaline[i,'genus_type'] %in% unique(animal$genus_type) || df.NonSaline[i,'genus_type'] %in% unique(saline$genus_type))
        df.NonSaline[i,'status'] = "Migrant"
      else
        df.NonSaline[i,'status'] = "Generalist"
}

df.NonSaline$status=as.factor(df.NonSaline$status)

################
############ Z-transform predictor
pvar='nb_genus'
datsc1 <- df.NonSaline
datsc1[pvar] <- lapply(df.NonSaline[pvar],scale)

datsc_nonsalin=datsc1
save(datsc_nonsalin,file="datsc_nonsalin.RData")


##full model 
glmer.nonsalin.1 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3/PI)+(nb_genus|sample),datsc1,family=poisson(link=log),
                            control=glmerControl(optimizer="bobyqa"))
summary(glmer.nonsalin.1)
save(glmer.nonsalin.1,file="glmm_nonS1.RData")


##supprimer la singularité
glmer.nonsalin.2 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|empo_3)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                            control=glmerControl(optimizer="bobyqa"))
save(glmer.nonsalin.2,file="glmm_nonsaline2.RData")
summary(glmer.nonsalin.2)

##supprimer empo_3 car variance = 0 (le plus significatif)
glmer.nonsalin.3 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                            control=glmerControl(optimizer="bobyqa"))
summary(glmer.nonsalin.3)
save(glmer.nonsalin.3,file="glmm_nonsaline3.RData")

##tester la significance des random d'après ben bolker : https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
#glmer.nonsalin.3 = glmer(nb_ASV~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
  #                          control=glmerControl(optimizer="bobyqa"))

glm1=update(glmer.nonsalin.3,.~nb_genus*status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI))
glm2=update(glmer.nonsalin.3,.~nb_genus*status+(nb_genus|genus_type/empo_3))
glm3 <- lm(nb_ASV~nb_genus*status,datsc1)
anova(glmer.nonsalin.3 ,glm1,glm2,glm3,glm4)

# Df    AIC    BIC logLik deviance     Chisq Chi Df Pr(>Chisq)    
# glm3             7 192209 192270 -96098   192195                                
# glm2            12 132101 132205 -66038   132077 60118.522      5  < 2.2e-16 ***
# glm1            15 131354 131484 -65662   131324   752.644      3  < 2.2e-16 ***
# glmm.nonsalin.3 16 131341 131479 -65654   131309    15.495      1  8.274e-05 ***

###significance of interaction
glmer.inter= glmer(nb_ASV~nb_genus+status+(nb_genus|genus_type/empo_3)+ (nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),
                                         control=glmerControl(optimizer="bobyqa"))

anova(glmm.nonsalin.3,glmer.inter)
# Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# glmer.inter     14 131408 131529 -65690   131380                             
# glmm.nonsalin.3 16 131341 131479 -65654   131309 71.222      2  3.422e-16 ***


#####
overdisp_fun(glmm.nonsalin.3)
# chisq        ratio          rdf            p 
# 2.501410e+04 5.946818e-01 4.206300e+04 1.000000e+00 
plot(fitted(glmm.nonsalin.3),residuals(glmm.nonsalin.3),xlab = "Fitted Values", ylab = "Residuals")
qqnorm(scale(resid(glmm.nonsalin.3)),ylab="Residual quantiles",col="orange")

##interacton plot
require(sjPlot)
pdf("Nonsaline.pdf")
plot1=plot_model(glmm.nonsalin.3,type="int",terms=c("nb_genus","status"),title="Non saline",show.legend = TRUE, axis.title=c("","ASV number/genus"),colors = c("black","red"))
plot1
dev.off()
