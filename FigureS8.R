
#### plots pour l'approche %identity
#### Figure S8

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
            axis.title = element_text(face = "bold",size = rel(0.8)),
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

##~/tableau_identity_1.0.tsv
tab100=read.table("~/tableau_identity_1.0.tsv",header=T,sep=",")
tab97=read.table("~/tableau_identity_0.97.tsv",header=T,sep=",")
tab95=read.table("~/tableau_identity_0.95.tsv",header=T,sep=",")
tab90=read.table("~/tableau_identity_0.9.tsv",header=T,sep=",")
tab85=read.table("~/tableau_identity_0.85.tsv",header=T,sep=",")
tab80=read.table("~/tableau_identity_0.8.tsv",header=T,sep=",")
tab75=read.table("~/tableau_identity_0.75.tsv",header=T,sep=",")

########
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

tab=tab100_97
tab100_97$nb_100=tab$nb_100/tab$nb_97

##################
##### linear, quadrati and cubic models 

### 80%~75%
fit_0=lm(nb_80 ~ nb_75, data = tab80_75)
fit_02=lm(nb_80 ~ poly(nb_75,3), data = tab80_75)
r_10=signif(summary(fit_0)$adj.r.squared, 2)
r_20=signif(summary(fit_02)$adj.r.squared, 2)

fitquad0=lm(nb_80 ~ poly(nb_75,2), data = tab80_75)
rquad0=signif(summary(fitquad0)$adj.r.squared, 2)

p_80_75=ggplot(fit_0$model, aes_string(x = names(fit_0$model)[2], y = names(fit_0$model)[1])) +
  geom_point(size=0.5) +
  stat_smooth(method = "lm", se=TRUE,col = "dimgrey")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 2, raw=TRUE),colour="blue")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 3, raw=TRUE),colour="red") +
  labs(x="75% clusters",y="80% clusters")+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r_10," ",rquad0," ",r_20),color="red", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r_10," ",rquad0),color="blue", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate( geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r_10),color="dimgrey", size = 3, 
            vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()

### 85%~80%
fit=lm(nb_85 ~ nb_80, data = tab80_85)
fit2=lm(nb_85 ~ poly(nb_80,3), data = tab80_85)

r1=signif(summary(fit)$adj.r.squared, 2)
r2=signif(summary(fit2)$adj.r.squared, 2)

fitquad=lm(nb_85 ~ poly(nb_80,2), data = tab80_85)
rquad=signif(summary(fitquad)$adj.r.squared, 2)

p85_80=ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
  geom_point(size=0.5) +
  stat_smooth(method = "lm", se=TRUE,col = "dimgrey")+
  stat_smooth(method = "lm", se=TRUE,formula=y ~ poly(x, 2, raw=TRUE),col = "blue")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 3, raw=TRUE),colour="red") +
  labs(x="80% clusters",y="85% clusters")+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r1," ",rquad," ",r2),color="red", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r1," ",rquad),color="blue", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate( geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r1),color="dimgrey", size = 3, 
            vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()

################# 90~85
fit00=lm(nb_90~nb_85,tab85_90)
fit002=lm(nb_90 ~ poly(nb_85,3), data = tab85_90)
r100=signif(summary(fit00)$adj.r.squared, 2)
r200=signif(summary(fit002)$adj.r.squared, 2)

fitquad00=lm(nb_90 ~ poly(nb_85,2), data = tab85_90)
rquad00=signif(summary(fitquad00)$adj.r.squared, 2)

pp_90_85=ggplot(fit00$model, aes_string(x = names(fit00$model)[2], y = names(fit00$model)[1])) +
  geom_point(size=0.5) +
  stat_smooth(method = "lm", se=TRUE,col = "dimgrey")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 2, raw=TRUE),colour="blue") +
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 3, raw=TRUE),colour="red") +
  labs(x="85% clusters",y="90% clusters")+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r100," ",rquad00," ",r200),color="red", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r100," ",rquad00),color="blue", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate( geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r100),color="dimgrey", size = 3, 
            vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()


####### 95%~90%
fit000=lm(nb_95~nb_90,tab95_9)
fit0002=lm(nb_95 ~ poly(nb_90,3), data = tab95_9)

r1000=signif(summary(fit000)$adj.r.squared, 2)
r2000=signif(summary(fit0002)$adj.r.squared, 2)

fitquad000=lm(nb_95 ~ poly(nb_90,2), data = tab95_9)
rquad000=signif(summary(fitquad000)$adj.r.squared, 2)

p_95_90=ggplot(fit000$model, aes_string(x = names(fit000$model)[2], y = names(fit000$model)[1])) +
  geom_point(size=0.5) +
  stat_smooth(method = "lm", se=TRUE,col = "dimgrey")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 2, raw=TRUE),colour="blue") +

  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 3, raw=TRUE),colour="red") +
  labs(x="90% clusters",y="95% clusters")+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r1000," ",rquad000," ",r2000),color="red", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r1000," ",rquad000),color="blue", size = 3, 
           vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate( geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r1000),color="dimgrey", size = 3, 
            vjust=1.2,hjust=0,fontface=2,fill = "white",label.size = 0)+

  theme_Publication()


### 97%~95%
fit0000=lm(nb_97~nb_95,tab95_97)
fit00002=lm(nb_97 ~ poly(nb_95,3), data = tab95_97)
r10000=signif(summary(fit0000)$adj.r.squared, 2)
r20000=signif(summary(fit00002)$adj.r.squared, 2)

fitquad0000=lm(nb_97 ~ poly(nb_95,2), data = tab95_97)
rquad0000=signif(summary(fitquad0000)$adj.r.squared, 2)

p_97_95=ggplot(fit0000$model, aes_string(x = names(fit0000$model)[2], y = names(fit0000$model)[1])) +
  geom_point(size=0.5) +
  stat_smooth(method = "lm", se=TRUE,col = "dimgrey")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 2, raw=TRUE),colour="blue") +

  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 3, raw=TRUE),colour="red") +
  labs(x="95% clusters",y="97% clusters")+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r10000," ",rquad0000," ",r20000),color="red", size = 3, 
           vjust=3,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r10000," ",rquad0000),color="blue", size = 3, 
           vjust=3,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate( geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r10000),color="dimgrey", size = 3, 
            vjust=3,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()


################ 100%~97%
fit_50=lm(nb_100~nb_97,tab100_97)
fit_502=lm(nb_100 ~ poly(nb_97,3), data = tab100_97)

r150=signif(summary(fit_50)$adj.r.squared, 2)
r152=signif(summary(fit_502)$adj.r.squared, 2)

fitquad50=lm(nb_100 ~ poly(nb_97,2), data = tab100_97)
rquad50=signif(summary(fitquad50)$adj.r.squared, 2)

p100_97=ggplot(fit_50$model, aes_string(x = names(fit_50$model)[2], y = names(fit_50$model)[1])) +
  geom_point(size=0.5) +
  stat_smooth(method = "lm", se=TRUE,col = "dimgrey")+
  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 2, raw=TRUE),colour="blue") +

  stat_smooth(method="lm",se=TRUE, formula=y ~ poly(x, 3, raw=TRUE),colour="red") +
  labs(x="97% clusters",y="100% clusters")+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r150," ",rquad50," ",r152),color="red", size = 3, 
           vjust=3,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate(geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r150," ",rquad50),color="blue", size = 3, 
           vjust=3,hjust=0,fontface=2,fill = "white",label.size = 0)+
  annotate( geom="label",x = -Inf, y = Inf,label = paste("Adj R2 : ",r150),color="dimgrey", size = 3, 
            vjust=3,hjust=0,fontface=2,fill = "white",label.size = 0)+
  theme_Publication()


#### plot
require(ggpubr)
pdf("FigureS8.pdf",width=8,heigh=6)
ggarrange(p_80_75,p85_80,pp_90_85, p_95_90,p_97_95,p100_97,
          #labels = c("A", "B","C","D","E","F"),
          ncol = 3, nrow = 2,
          font.label = list(size = 11))
dev.off()
