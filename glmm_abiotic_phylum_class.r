#### Abiotic factors analysis
#### glmm Class:Phylum
#### NM

library(lme4)
library("readxl")
library(ggplot2)
library(sjPlot)

df.metadata=read.csv("/Users/naima/Projet_Diversification/Metadata/emp_qiime_mapping_subset_2k.tsv",sep="\t")

df.subset=df.metadata[,c("X.SampleID","latitude_deg", "elevation_m","temperature_deg_c", "ph","empo_3")]

hist(df.subset$ph)
hist(df.subset$latitude_deg)
hist(df.subset$elevation_m)
hist(df.subset$temperature_deg_c)

subset.temp.ph=df.subset[which(!is.na(df.subset$temperature_deg_c) & !is.na(df.subset$ph) & !is.na(df.subset$latitude_deg) & !is.na(df.subset$elevation_m)),]
dim(subset.temp.ph)
head(subset.temp.ph)
length(unique(subset.temp.ph$X.SampleID))
dim(subset.temp.ph[which(is.na(subset.temp.ph$latitude_deg)),])

hist(subset.temp.ph$ph)
hist(subset.temp.ph$latitude_deg)
hist(subset.temp.ph$elevation_m)
hist(subset.temp.ph$temperature_deg_c)

### this table contains phyla count (diversity) and class:phyla ratio (diversification) for every phylum in all samples
### could be a diversity-diversification table for any other ratio.

tab_phyla_class=read.table("/Users/naima/Projet_Diversification/Diversification_rate_18mars/phyla_tab_2k.tsv",header=T,sep=",")
head(tab_phyla_class)
df.sample=unique(tab_phyla_class$sample)
length(df.sample)

study_lab=read_excel("/Users/naima/Projet_Diversification/Diversification_rate_11mars/lab_emp.xlsx", col_names  = TRUE)
###lire les environnement emp1,emp2,emp3
dim(study_lab)
colnames(study_lab)=c("study","PI")

for (i in 1:nrow(tab_phyla_class) )
  tab_phyla_class[i,'study']=strsplit(as.character(tab_phyla_class[i,]$sample),"[.]")[[1]][1]

df.study=merge(tab_phyla_class,study_lab,by="study")

colnames(subset.temp.ph)[1]="sample"

df.data=merge(df.study,subset.temp.ph,by="sample")

colnames(df.data)
head(df.data)

df.data$empo_3=as.factor(df.data$empo_3)
df.data$phylum_type=as.factor(df.data$phylum_type)
df.data$PI=as.factor(df.data$PI)

range(df.data$ph)
range(df.data$latitude_deg)
range(df.data$elevation_m)
range(df.data$temperature_deg_c)

### Z-transformation
pvar='nb_phylum'
pvar1='ph'
pvar2="temperature_deg_c"
pvar3="latitude_deg"
pvar4="elevation_m"

datsc <- df.data
datsc[pvar] <- lapply(df.data[pvar],scale)
datsc[pvar1] <- lapply(df.data[pvar1],scale)
datsc[pvar2] <- lapply(df.data[pvar2],scale)
datsc[pvar3] <- lapply(df.data[pvar3],scale)
datsc[pvar4] <- lapply(df.data[pvar4],scale)


glmer1.phylum.intercept=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3/PI)+(1|sample),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
summary(glmer1.phylum.intercept)

##enlever la singularité sur empo3 et sample
g1=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3:PI),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
summary(g1)

## tester le segond degrés des facteurs abiotiques
g1.temp=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+I(temperature_deg_c^2)+nb_phylum:I(temperature_deg_c^2)+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3:PI),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(g1.temp,g1) ##poly signif mais only 4 de diff entre AIC alors on garde g1

g1.elev=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+I(elevation_m^2)+nb_phylum:I(elevation_m^2)+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3:PI),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(g1.elev,g1) 

g1.lat=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+I(latitude_deg^2)+nb_phylum:I(latitude_deg^2)+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3:PI),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(g1.lat,g1) 

g1.ph=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+I(ph^2)+nb_phylum:I(ph^2)+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3:PI),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(g1.ph,g1) 

### g1 is significant 
save(g1,file="g1_signif.RData")

## Random Effect significance
g1=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3)+(1|empo_3:PI),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
g2=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
g3=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
g4=lm(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m,datsc)
anova(g1,g2,g3,g4)  ##g2 est significant

summary(g2)
save(g2,file="g2_sign.RData")

## fixed effect significance : anova between full and null model (without fixed effect)
#g2=glmer(nb_class~nb_phylum+temperature_deg_c+latitude_deg+ph+elevation_m+nb_phylum:temperature_deg_c+nb_phylum:latitude_deg+nb_phylum:ph+nb_phylum:elevation_m+(1|phylum_type/empo_3),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
null_m=glmer(nb_class~1+(1|phylum_type/empo_3),datsc,family=poisson(link=log),control=glmerControl(optimizer="bobyqa"))
anova(g2,null_m)
