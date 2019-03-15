##### GLMM  Class:Phylum ratio
##### 28 oct NM

if(!require(lme4)){install.packages("lme4")}
require(lme4)

if(!require(readxl)){install.packages("readxl")}
require(readxl)

if(!require(lattice)){install.packages("lattice")}
require(lattice)

if(!require(jtools)){install.packages("jtools")}
require(jtools)

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
colnames(tab_emp)
colnames(tab_emp)=c("sample","empo_0", "empo_1" ,"empo_2","empo_3")

tab_emp=as.data.frame(tab_emp)
tab_emp$empo_1=as.factor(tab_emp$empo_1)
tab_emp$empo_2=as.factor(tab_emp$empo_2)
tab_emp$empo_3=as.factor(tab_emp$empo_3)

###table des nb_phyla, phyla type et le nombre correspondant en classe pour chaque sample
tab_phylatype=read.table("~/phyla_tab_2k.tsv",header=T,sep=",")
tab_phylatype=as.data.frame(tab_phylatype)
tab_phylatype$phylum_type=as.factor(tab_phylatype$phylum_type)
tab_phylatype_emp=merge(tab_phylatype[,-1],tab_emp[,c("sample","empo_3")],by="sample")

### ajouter le lab
for (i in 1:nrow(tab_phylatype_emp) )
  tab_phylatype_emp[i,'study']=strsplit(as.character(tab_phylatype_emp[i,]$sample),"[.]")[[1]][1]

tab.final=merge(tab_phylatype_emp,study_lab,by="study")
tab.final$PI=as.factor(tab.final$PI)
tab.final$phylum_type=as.factor(tab.final$phylum_type)
tab.final$empo_3=as.factor(tab.final$empo_3)

######################
save(tab.final,file="tab_final_phylum.RData")
##########

hist(tab.final$nb_class)
hist(tab.final$nb_phylum)

####   z transform the predictor
pvar='nb_phylum'
datsc <- tab.final
datsc[pvar] <- lapply(tab.final[pvar],scale)

##full model
glmer.phylum = glmer(nb_class~nb_phylum+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3/PI)+(nb_phylum|sample),datsc,family=poisson(link=log),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
save(glmer.phylum,file="glmer_phylum.RData")
summary(glmer.phylum)

##### supprimer la singularité : sample variation is 0
glmer.phylum.sample.1 = glmer(nb_class~nb_phylum+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3/PI),datsc,family=poisson(link=log),
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

save(glmer.phylum.sample.1,file="glmer_phylum_sample_sept.RData")
summary(glmer.phylum.sample.1)

##### supprimer la singularité: il y a une correlation de 1 entre slope et intercept de biome, on garde biome juste sur slope
glmer.phylum.biomSlop = glmer(nb_class~nb_phylum+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3:PI)+(nb_phylum-1|empo_3),datsc,family=poisson(link=log),
                              control=glmerControl(optimizer="bobyqa"))
summary(glmer.phylum.biomSlop)
save(glmer.phylum.biomSlop,file="glmer_phylum_biomSlop.RData")

###tester la significance des randoms
glmer.phylum.biomSlop1 = glmer(nb_class~nb_phylum+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3:PI),datsc,family=poisson(link=log),
                              control=glmerControl(optimizer="bobyqa"))

glmer.phylum.biomSlop2 = glmer(nb_class~nb_phylum+(nb_phylum|phylum_type/empo_3),datsc,family=poisson(link=log),
                              control=glmerControl(optimizer="bobyqa"))

glmer.phylum.biomSlop3 = lm(nb_class~nb_phylum,datsc)

anova(glmer.phylum.sample.1,glmer.phylum.biomSlop,glmer.phylum.biomSlop1,glmer.phylum.biomSlop2,glmer.phylum.biomSlop3)
##### glmer.phylum.biomSlop est le plus significatif

#### test the fixed effect significance
g1=update(glmer.phylum.biomSlop,.~1+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3:PI)+(nb_phylum-1|empo_3))
anova(glmer.phylum.biomSlop,g1)  ## glmer.phylum.biomSlop is significant

##### overdispersion test and diagnostic plots
overdisp_fun(glmer.phylum.biomSlop)
plot(glmer.phylum.biomSlop)
qqnorm(residuals(glmer.phylum.biomSlop))

######
#### DBD variation across environments
#### biome as fixed effect
#### model without the main effects (of empo3 and nb_phylum, i.e. separate intercept
#### for every biome)
##########
glmer.phylum.biom.2 = glmer(nb_class~nb_phylum:empo_3+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3:PI),datsc,family=poisson(link=log),
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
summary(glmer.phylum.biom.2)
save(glmer.phylum.biom.2,file="glmer_phylum_biomFixed2.RData")

#####significance de l'interaction
g_1 = glmer(nb_class~1+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3:PI),datsc,family=poisson(link=log),
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
anova(glmer.phylum.biom.2,g_1)

##### interaction + the main effects of nb_phylum and empo_3
glmer.phylum.biom.21 = glmer(nb_class~nb_phylum*empo_3+(nb_phylum|phylum_type/empo_3)+(nb_phylum|empo_3:PI),datsc,family=poisson(link=log),
                            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
save(glmer.phylum.biom.21,file="glmerphylumbiom21.RData")  ##only 4 aic difference
anova(glmer.phylum.biom.2,glmer.phylum.biom.21) ## dAIC = 4, weak dAIC means the two models are equivalent, we keep the parcimonious model
##

#### diagnostic plots and overdispersion
interact_plot(glmer.phylum.biom.2,pred=nb_phylum,modx=empo_3)
overdisp_fun(glmer.phylum.biom.2)
plot(glmer.phylum.biom.2)
qqnorm(residuals(glmer.phylum.biom.2))

