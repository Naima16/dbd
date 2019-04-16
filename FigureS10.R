
#### plot 85%~80% by biome (rsquared of lm deg 1,2,3)
#### 28 janv 2018

require(plyr)
require(ggplot2)

theme_Publication <- function(base_size=14, base_family="Helvetica"){ #Comic Sans MS"){ ##helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size=8), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold",size = 11)
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
tab100=read.table("~/tableau_identity_1.0.tsv",header=T,sep=",")
tab97=read.table("~/tableau_identity_0.97.tsv",header=T,sep=",")
tab95=read.table("~/tableau_identity_0.95.tsv",header=T,sep=",")
tab90=read.table("~/tableau_identity_0.9.tsv",header=T,sep=",")
tab85=read.table("~/tableau_identity_0.85.tsv",header=T,sep=",")
tab80=read.table("~/tableau_identity_0.8.tsv",header=T,sep=",")
tab75=read.table("~/tableau_identity_0.75.tsv",header=T,sep=",")

colnames(tab75)[3]='nb_75'
colnames(tab80)[3]='nb_80'
colnames(tab85)[3]='nb_85'
colnames(tab90)[3]='nb_90'
colnames(tab95)[3]='nb_95'
colnames(tab97)[3]='nb_97'
colnames(tab100)[3]='nb_100'

tab80_75=merge(tab80[,-1],tab75[,-1],by='sample')
tab80_85=merge(tab80[,-1],tab85[,-1],by='sample')
tab85_90=merge(tab85[,-1],tab90[,-1],by='sample')
tab95_9=merge(tab95[,-1],tab90[,-1],by='sample')
tab95_97=merge(tab95[,-1],tab97[,-1],by='sample')
tab100_97=merge(tab97[,-1],tab100[,-1],by='sample')

tab=tab100_97
tab100_97$nb_100=tab$nb_100/tab$nb_97

tab=tab80_75
tab80_75$nb_80=tab$nb_80/tab$nb_75

tab=tab80_85
tab80_85$nb_85=tab$nb_85/tab$nb_80

tab=tab85_90
tab85_90$nb_90=tab$nb_90/tab$nb_85

tab=tab95_9
tab95_9$nb_95=tab$nb_95/tab$nb_90

head(tab95_97)
tab=tab95_97
tab95_97$nb_97=tab$nb_97/tab$nb_95


#### biome assignation to samples
tab_emp=read.csv("~/sample_emplevel.tsv",header=T)
head(tab_emp)
tab_emp=tab_emp[,-1]
colnames(tab_emp)=c("sample","empo_0", "empo_1" ,    "empo_2"  ,   "empo_3")

tab_emp=as.data.frame(tab_emp)
tab_emp$empo_1=as.factor(tab_emp$empo_1)
tab_emp$empo_2=as.factor(tab_emp$empo_2)
tab_emp$empo_3=as.factor(tab_emp$empo_3)


###
tab80_75.emp=merge(tab80_75,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")
tab85_80.emp=merge(tab80_85,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")
tab90_85.emp=merge(tab85_90,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")
tab95_90.emp=merge(tab95_9,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")
tab95_97.emp=merge(tab95_97,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")
tab97_100.emp=merge(tab100_97,tab_emp[,c("sample","empo_1","empo_2","empo_3")],by="sample")


###

################################## plot les modèles lineaire, quadratic et cubic
###80~75
mydata=tab80_75.emp

mods <- dlply(mydata, .(empo_3), function(x) lm(nb_80 ~ nb_75, data = x))
mods_2 <- dlply(mydata, .(empo_3), function(x) lm(nb_80 ~ poly(nb_75,2), data = x))
mods_3 <- dlply(mydata, .(empo_3), function(x) lm(nb_80 ~ poly(nb_75,3), data = x))

# Functions to extract the slope and p-value from a generic lm model object x

rsq=function(x) summary(x)$adj.r.squared
pval <- function(x) p.adjust(anova(x)$`Pr(>F)`[1],"bonferroni")

# Apply the 3 functions to the list of models - both yield two column
# data frames with individual as one of the variables. This is convenient
pval1=ldply(mods, pval) 
colnames(pval1)=c("empo_3","pv1")
pval2=ldply(mods_2, pval) 
colnames(pval2)=c("empo_3","pv2")
pval3=ldply(mods_3, pval) 
colnames(pval3)=c("empo_3","pv3")

s3=ldply(mods, rsq)
colnames(s3)=c("empo_3","r1")
s4=ldply(mods_2, rsq)
colnames(s4)=c("empo_3","r2")
s5=ldply(mods_3, rsq)
colnames(s5)=c("empo_3","r3")

### merge all together 
ds_out=Reduce(function(...) merge(..., all = TRUE, by = "empo_3"),
              list(s3,s4,s5,pval1,pval2,pval3))

ds_out$rlab1 <- round(ds_out$r1, 2) 
ds_out$rlab2 <- round(ds_out$r2, 2)
ds_out$rlab3 <- round(ds_out$r3, 2)

ds_out$sig1 <- ifelse(ds_out$pv1 <= 0.05, 1, 0)
ds_out$sig2 <- ifelse(ds_out$pv2 <= 0.05, 1, 0)
ds_out$sig3 <- ifelse(ds_out$pv3 <= 0.05, 1, 0)


## plot les 3 modèles lineaire, cuadr et cub
dd=merge(mydata,ds_out[,c("empo_3","rlab1","rlab2","rlab3","sig1","sig2","sig3")],by="empo_3")

pdf("plot_75_80.pdf",width=10.2,height=10.2)

Plot <- ggplot(data=dd,aes(nb_75,nb_80),color = factor(empo_3))
Plot_ByBiom<-Plot + geom_point(size=1) + xlab("75% clusters") + ylab(paste("80% clusters","")) +
  geom_smooth(aes(linetype=factor(sig1)),method="lm",se=FALSE,fullrange=FALSE,color="dimgrey", show.legend = F)+
  stat_smooth(aes(linetype=factor(sig2)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 2, raw=TRUE),colour="blue", show.legend = F) +
  stat_smooth(aes(linetype=factor(sig3)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 3, raw=TRUE),colour="red", show.legend = F) +
  scale_linetype_manual(values=c("twodash","solid"))+
  facet_wrap( ~ empo_3, scales = "free")+
  geom_label(data = ds_out, aes(x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2," ",rlab3)),show.legend=F,color="red", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p1=Plot_ByBiom + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2)),show.legend=F,color="blue", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p80_75= p1 + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1)),show.legend=F,color="dimgrey", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()
p80_75
dev.off()

###85~80
 mydata=tab85_80.emp
 head(mydata)
 range(mydata$nb_85)
 range(mydata$nb_80)

mods <- dlply(mydata, .(empo_3), function(x) lm(nb_85 ~ nb_80, data = x))
mods_2 <- dlply(mydata, .(empo_3), function(x) lm(nb_85 ~ poly(nb_80,2), data = x))
mods_3 <- dlply(mydata, .(empo_3), function(x) lm(nb_85 ~ poly(nb_80,3), data = x))
 
# Functions to extract the slope and p-value from a generic lm model object x

rsq=function(x) summary(x)$adj.r.squared
pval <- function(x) p.adjust(anova(x)$`Pr(>F)`[1],"bonferroni")
 
# Apply the 3 functions to the list of models - both yield two column
# data frames with individual as one of the variables. This is convenient
pval1=ldply(mods, pval) 
colnames(pval1)=c("empo_3","pv1")
pval2=ldply(mods_2, pval) 
colnames(pval2)=c("empo_3","pv2")
pval3=ldply(mods_3, pval) 
colnames(pval3)=c("empo_3","pv3")

s3=ldply(mods, rsq)
colnames(s3)=c("empo_3","r1")
s4=ldply(mods_2, rsq)
colnames(s4)=c("empo_3","r2")
s5=ldply(mods_3, rsq)
colnames(s5)=c("empo_3","r3")

### merge all together 
ds_out=Reduce(function(...) merge(..., all = TRUE, by = "empo_3"),
       list(s3,s4,s5,pval1,pval2,pval3))

ds_out$rlab1 <- round(ds_out$r1, 2) 
ds_out$rlab2 <- round(ds_out$r2, 2)
ds_out$rlab3 <- round(ds_out$r3, 2)

ds_out$sig1 <- ifelse(ds_out$pv1 <= 0.05, 1, 0)
ds_out$sig2 <- ifelse(ds_out$pv2 <= 0.05, 1, 0)
ds_out$sig3 <- ifelse(ds_out$pv3 <= 0.05, 1, 0)


## plot les 3 modèles lineaire, cuadr et cub
dd=merge(mydata,ds_out[,c("empo_3","rlab1","rlab2","rlab3","sig1","sig2","sig3")],by="empo_3")

pdf("plot_85_80.pdf",width=10.2,height=10.2)

Plot <- ggplot(data=dd,aes(nb_80,nb_85),color = factor(empo_3))
Plot_ByBiom<-Plot + geom_point(size=1) + xlab("80% clusters") + ylab(paste("85% clusters","")) +
  geom_smooth(aes(linetype=factor(sig1)),method="lm",se=FALSE,fullrange=FALSE,color="dimgrey", show.legend = F)+
  stat_smooth(aes(linetype=factor(sig2)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 2, raw=TRUE),colour="blue", show.legend = F) +
  stat_smooth(aes(linetype=factor(sig3)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 3, raw=TRUE),colour="red", show.legend = F) +
  scale_linetype_manual(values=c("twodash","solid"))+
  facet_wrap( ~ empo_3, scales = "free")+
  geom_label(data = ds_out, aes(x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2," ",rlab3)),show.legend=F,color="red", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p1=Plot_ByBiom + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2)),show.legend=F,color="blue", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p85_80= p1 + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1)),show.legend=F,color="dimgrey", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
    theme_Publication()
p85_80
dev.off()

### 90~85
mydata=tab90_85.emp

mods <- dlply(mydata, .(empo_3), function(x) lm(nb_90 ~ nb_85, data = x))
mods_2 <- dlply(mydata, .(empo_3), function(x) lm(nb_90 ~ poly(nb_85,2), data = x))
mods_3 <- dlply(mydata, .(empo_3), function(x) lm(nb_90 ~ poly(nb_85,3), data = x))

# Functions to extract the slope and p-value from a generic lm model object x

rsq=function(x) summary(x)$adj.r.squared
pval <- function(x) p.adjust(anova(x)$`Pr(>F)`[1],"bonferroni")

# Apply the 3 functions to the list of models - both yield two column
# data frames with individual as one of the variables. This is convenient
pval1=ldply(mods, pval) 
colnames(pval1)=c("empo_3","pv1")
pval2=ldply(mods_2, pval) 
colnames(pval2)=c("empo_3","pv2")
pval3=ldply(mods_3, pval) 
colnames(pval3)=c("empo_3","pv3")

s3=ldply(mods, rsq)
colnames(s3)=c("empo_3","r1")
s4=ldply(mods_2, rsq)
colnames(s4)=c("empo_3","r2")
s5=ldply(mods_3, rsq)
colnames(s5)=c("empo_3","r3")

### merge all together 
ds_out=Reduce(function(...) merge(..., all = TRUE, by = "empo_3"),
              list(s3,s4,s5,pval1,pval2,pval3))

ds_out$rlab1 <- round(ds_out$r1, 2) 
ds_out$rlab2 <- round(ds_out$r2, 2)
ds_out$rlab3 <- round(ds_out$r3, 2)

ds_out$sig1 <- ifelse(ds_out$pv1 <= 0.05, 1, 0)
ds_out$sig2 <- ifelse(ds_out$pv2 <= 0.05, 1, 0)
ds_out$sig3 <- ifelse(ds_out$pv3 <= 0.05, 1, 0)

## plot les 3 modèles lineaire, cuadr et cub
dd=merge(mydata,ds_out[,c("empo_3","rlab1","rlab2","rlab3","sig1","sig2","sig3")],by="empo_3")

pdf("plot_90_85.pdf",width=10.2,height=10.2)

Plot <- ggplot(data=dd,aes(nb_85,nb_90),color = factor(empo_3))
Plot_ByBiom<-Plot + geom_point(size=1) + xlab("85% clusters") + ylab(paste("90% clusters","")) +
  geom_smooth(aes(linetype=factor(sig1)),method="lm",se=FALSE,fullrange=FALSE,color="dimgrey", show.legend = F)+
  stat_smooth(aes(linetype=factor(sig2)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 2, raw=TRUE),colour="blue", show.legend = F) +
  stat_smooth(aes(linetype=factor(sig3)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 3, raw=TRUE),colour="red", show.legend = F) +
  scale_linetype_manual(values=c("twodash","solid"))+
  facet_wrap( ~ empo_3, scales = "free")+
  geom_label(data = ds_out, aes(x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2," ",rlab3)),show.legend=F,color="red", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p1=Plot_ByBiom + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2)),show.legend=F,color="blue", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p90_85= p1 + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1)),show.legend=F,color="dimgrey", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()
p90_85
dev.off()

##### 95~90
mydata=tab95_90.emp

mods <- dlply(mydata, .(empo_3), function(x) lm(nb_95 ~ nb_90, data = x))
mods_2 <- dlply(mydata, .(empo_3), function(x) lm(nb_95 ~ poly(nb_90,2), data = x))
mods_3 <- dlply(mydata, .(empo_3), function(x) lm(nb_95 ~ poly(nb_90,3), data = x))

# Functions to extract the slope and p-value from a generic lm model object x

rsq=function(x) summary(x)$adj.r.squared
pval <- function(x) p.adjust(anova(x)$`Pr(>F)`[1],"bonferroni")

# Apply the 3 functions to the list of models - both yield two column
# data frames with individual as one of the variables. This is convenient
pval1=ldply(mods, pval) 
colnames(pval1)=c("empo_3","pv1")
pval2=ldply(mods_2, pval) 
colnames(pval2)=c("empo_3","pv2")
pval3=ldply(mods_3, pval) 
colnames(pval3)=c("empo_3","pv3")

s3=ldply(mods, rsq)
colnames(s3)=c("empo_3","r1")
s4=ldply(mods_2, rsq)
colnames(s4)=c("empo_3","r2")
s5=ldply(mods_3, rsq)
colnames(s5)=c("empo_3","r3")

### merge all together 
ds_out=Reduce(function(...) merge(..., all = TRUE, by = "empo_3"),
              list(s3,s4,s5,pval1,pval2,pval3))

ds_out$rlab1 <- round(ds_out$r1, 2) 
ds_out$rlab2 <- round(ds_out$r2, 2)
ds_out$rlab3 <- round(ds_out$r3, 2)

ds_out$sig1 <- ifelse(ds_out$pv1 <= 0.05, 1, 0)
ds_out$sig2 <- ifelse(ds_out$pv2 <= 0.05, 1, 0)
ds_out$sig3 <- ifelse(ds_out$pv3 <= 0.05, 1, 0)

## plot les 3 modèles lineaire, cuadr et cub
dd=merge(mydata,ds_out[,c("empo_3","rlab1","rlab2","rlab3","sig1","sig2","sig3")],by="empo_3")

pdf("plot_95_90.pdf",width=10.2,height=10.2)

Plot <- ggplot(data=dd,aes(nb_90,nb_95),color = factor(empo_3))
Plot_ByBiom<-Plot + geom_point(size=1) + xlab("90% clusters") + ylab(paste("95% clusters","")) +
  geom_smooth(aes(linetype=factor(sig1)),method="lm",se=FALSE,fullrange=FALSE,color="dimgrey", show.legend = F)+
  stat_smooth(aes(linetype=factor(sig2)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 2, raw=TRUE),colour="blue", show.legend = F) +
  stat_smooth(aes(linetype=factor(sig3)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 3, raw=TRUE),colour="red", show.legend = F) +
  scale_linetype_manual(values=c("twodash","solid"))+
  facet_wrap( ~ empo_3, scales = "free")+
  geom_label(data = ds_out, aes(x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2," ",rlab3)),show.legend=F,color="red", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p1=Plot_ByBiom + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2)),show.legend=F,color="blue", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p95_90= p1 + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1)),show.legend=F,color="dimgrey", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()
p95_90
dev.off()

##### 97~95
mydata=tab95_97.emp

mods <- dlply(mydata, .(empo_3), function(x) lm(nb_97 ~ nb_95, data = x))
mods_2 <- dlply(mydata, .(empo_3), function(x) lm(nb_97 ~ poly(nb_95,2), data = x))
mods_3 <- dlply(mydata, .(empo_3), function(x) lm(nb_97 ~ poly(nb_95,3), data = x))

# Functions to extract the slope and p-value from a generic lm model object x

rsq=function(x) summary(x)$adj.r.squared
pval <- function(x) p.adjust(anova(x)$`Pr(>F)`[1],"bonferroni")

# Apply the 3 functions to the list of models - both yield two column
# data frames with individual as one of the variables. This is convenient
pval1=ldply(mods, pval) 
colnames(pval1)=c("empo_3","pv1")
pval2=ldply(mods_2, pval) 
colnames(pval2)=c("empo_3","pv2")
pval3=ldply(mods_3, pval) 
colnames(pval3)=c("empo_3","pv3")

s3=ldply(mods, rsq)
colnames(s3)=c("empo_3","r1")
s4=ldply(mods_2, rsq)
colnames(s4)=c("empo_3","r2")
s5=ldply(mods_3, rsq)
colnames(s5)=c("empo_3","r3")

### merge all together 
ds_out=Reduce(function(...) merge(..., all = TRUE, by = "empo_3"),
              list(s3,s4,s5,pval1,pval2,pval3))

ds_out$rlab1 <- round(ds_out$r1, 2) 
ds_out$rlab2 <- round(ds_out$r2, 2)
ds_out$rlab3 <- round(ds_out$r3, 2)

ds_out$sig1 <- ifelse(ds_out$pv1 <= 0.05, 1, 0)
ds_out$sig2 <- ifelse(ds_out$pv2 <= 0.05, 1, 0)
ds_out$sig3 <- ifelse(ds_out$pv3 <= 0.05, 1, 0)

## plot les 3 modèles lineaire, cuadr et cub
dd=merge(mydata,ds_out[,c("empo_3","rlab1","rlab2","rlab3","sig1","sig2","sig3")],by="empo_3")

pdf("plot_97_95.pdf",width=10.2,height=10.2)

Plot <- ggplot(data=dd,aes(nb_95,nb_97),color = factor(empo_3))
Plot_ByBiom<-Plot + geom_point(size=1) + xlab("95% clusters") + ylab(paste("97% clusters","")) +
  geom_smooth(aes(linetype=factor(sig1)),method="lm",se=FALSE,fullrange=FALSE,color="dimgrey", show.legend = F)+
  stat_smooth(aes(linetype=factor(sig2)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 2, raw=TRUE),colour="blue", show.legend = F) +
  stat_smooth(aes(linetype=factor(sig3)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 3, raw=TRUE),colour="red", show.legend = F) +
  scale_linetype_manual(values=c("twodash","solid"))+
  facet_wrap( ~ empo_3, scales = "free")+
  geom_label(data = ds_out, aes(x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2," ",rlab3)),show.legend=F,color="red", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p1=Plot_ByBiom + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2)),show.legend=F,color="blue", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p97_95= p1 + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1)),show.legend=F,color="dimgrey", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()
p97_95
dev.off()

### 100~97
mydata=tab97_100.emp

mods <- dlply(mydata, .(empo_3), function(x) lm(nb_100 ~ nb_97, data = x))
mods_2 <- dlply(mydata, .(empo_3), function(x) lm(nb_100 ~ poly(nb_97,2), data = x))
mods_3 <- dlply(mydata, .(empo_3), function(x) lm(nb_100 ~ poly(nb_97,3), data = x))

# Functions to extract the slope and p-value from a generic lm model object x

rsq=function(x) summary(x)$adj.r.squared
pval <- function(x) p.adjust(anova(x)$`Pr(>F)`[1],"bonferroni")

# Apply the 3 functions to the list of models - both yield two column
# data frames with individual as one of the variables. This is convenient
pval1=ldply(mods, pval) 
colnames(pval1)=c("empo_3","pv1")
pval2=ldply(mods_2, pval) 
colnames(pval2)=c("empo_3","pv2")
pval3=ldply(mods_3, pval) 
colnames(pval3)=c("empo_3","pv3")

s3=ldply(mods, rsq)
colnames(s3)=c("empo_3","r1")
s4=ldply(mods_2, rsq)
colnames(s4)=c("empo_3","r2")
s5=ldply(mods_3, rsq)
colnames(s5)=c("empo_3","r3")

### merge all together 
ds_out=Reduce(function(...) merge(..., all = TRUE, by = "empo_3"),
              list(s3,s4,s5,pval1,pval2,pval3))

ds_out$rlab1 <- round(ds_out$r1, 2) 
ds_out$rlab2 <- round(ds_out$r2, 2)
ds_out$rlab3 <- round(ds_out$r3, 2)

ds_out$sig1 <- ifelse(ds_out$pv1 <= 0.05, 1, 0)
ds_out$sig2 <- ifelse(ds_out$pv2 <= 0.05, 1, 0)
ds_out$sig3 <- ifelse(ds_out$pv3 <= 0.05, 1, 0)

## plot les 3 modèles lineaire, cuadr et cub
dd=merge(mydata,ds_out[,c("empo_3","rlab1","rlab2","rlab3","sig1","sig2","sig3")],by="empo_3")

pdf("plot_97_100.pdf",width=10.2,height=10.2)

Plot <- ggplot(data=dd,aes(nb_97,nb_100),color = factor(empo_3))
Plot_ByBiom<-Plot + geom_point(size=1) + xlab("97% clusters") + ylab(paste("100% clusters","")) +
  geom_smooth(aes(linetype=factor(sig1)),method="lm",se=FALSE,fullrange=FALSE,color="dimgrey", show.legend = F)+
  stat_smooth(aes(linetype=factor(sig2)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 2, raw=TRUE),colour="blue", show.legend = F) +
  stat_smooth(aes(linetype=factor(sig3)),method="lm",se=FALSE, fullrange=FALSE,formula=y ~ poly(x, 3, raw=TRUE),colour="red", show.legend = F) +
 
  scale_linetype_manual(values=c("twodash","solid"))+
  facet_wrap( ~ empo_3, scales = "free")+
  geom_label(data = ds_out, aes(x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2," ",rlab3)),show.legend=F,color="red", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p1=Plot_ByBiom + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1," ",rlab2)),show.legend=F,color="blue", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)
p100_97= p1 + geom_label(data = ds_out, aes( x = -Inf, y = Inf,label = paste("Adj R2 : ",rlab1)),show.legend=F,color="dimgrey", size = 3, vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()
p100_97
dev.off()

