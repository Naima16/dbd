#### glmer genus ASV avec genome size en variable continue
#### glmer avec genome size

library(lme4)
library("readxl")
library(sjPlot)

### table des PI
study_lab=read_excel("~/lab_emp.xlsx", col_names  = TRUE)
###lire les environnement emp1,emp2,emp3
dim(study_lab)
colnames(study_lab)=c("study","PI")
head(study_lab)

##### table des empo pour chaque sample
tab_emp=read.csv("~/sample_emplevel.tsv",header=T)
dim(tab_emp)
tab.emp3=table(tab_emp$empo_3)
as.data.frame(tab.emp3)
tab_emp=tab_emp[,-1]
colnames(tab_emp)
colnames(tab_emp)=c("sample","empo_0", "empo_1" ,    "empo_2"  ,   "empo_3")
tab_emp[1,"sample"]
head(tab_emp)

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
########## genus genome size

genome.size=read_excel("~/genome_sizes_NM.xlsx")
gen.size=genome.size[,1:2]
colnames(gen.size)=c("genus_type","size")

data.gensize=merge(tab.final,gen.size,by="genus_type")

pvar='nb_genus'
pvar2='size'
datsc1 <- data.gensize

datsc1[pvar] <- lapply(data.gensize[pvar],scale)
datsc1[pvar2] <- lapply(data.gensize[pvar2],scale)
save(datsc1,file="datsc1.RData")

######### modèle glmer avec size

glmer.size.1 = glmer(nb_ASV~nb_genus*size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3/PI)+(nb_genus|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
summary(glmer.size.1)
save(glmer.size.1,file="glmer_size_1.RData")

## enlever la singularité empo3 et sample du model glmer.size.3
glmer.size.2 = glmer(nb_ASV~nb_genus*size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3:PI)+(nb_genus-1|empo_3)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
summary(glmer.size.2)
save(glmer.size.2,file="glmersize2.RData")

## enlever empo3 car var = 0
glmer.size.3 = glmer(nb_ASV~nb_genus*size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
summary(glmer.size.3)
save(glmer.size.3,file="glmersize3.RData")


###le modèle glmer.size.3 avec seulement interaction (sans main effects)
glmer.size.4 = glmer(nb_ASV~nb_genus:size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
anova(glmer.size.3,glmer.size.4) #3 est significant
save(glmer.size.4,file="glmersize4.RData")


## tester la significance de l'interaction
glmer.size.5 = glmer(nb_ASV~nb_genus+size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
anova(glmer.size.3,glmer.size.5)

#### tester les random effects
#glmer.size.3 = glmer(nb_ASV~nb_genus*size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3:PI)+(nb_genus-1|sample),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,

glm1 = glmer(nb_ASV~nb_genus*size+(nb_genus|genus_type/empo_3)+(nb_genus|empo_3:PI),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
glm2 = glmer(nb_ASV~nb_genus*size+(nb_genus|genus_type/empo_3),datsc1,family=poisson(link=log),control=glmerControl(optimizer="bobyqa")) #,
glm3 = lm(nb_ASV~nb_genus*size,data=datsc1)
anova(glmer.size.3 ,glm1,glm2,glm3)

overdisp_fun(glmer.size.3)

## diagnostic plots
plot(fitted(glmer.size.3),residuals(glmer.size.3),xlab = "Fitted Values", ylab = "Residuals")
qqnorm(scale(resid(glmer.size.3)),ylab="Residual quantiles",col="orange")


## plot interaction
johnson_neyman(glmer.size.3,pred=nb_genus,modx=size)
interact_plot(glmer.size.3,pred=nb_genus,modx=size,interval = T)
