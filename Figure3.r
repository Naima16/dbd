
#### Figure 3 , PCA and indicator species analysis
#### NM 8 avril 2019


if(!require(mvpart)){
  library(devtools)
  install_github("cran/mvpart", force = TRUE)}

if(!require(lme4)){install.packages("lme4")}
require(lme4)

if(!require(readxl)){install.packages("readxl")}
require(readxl)

if(!require(vegan)){install.packages("vegan")} ##rda and decostand
require(vegan)

if(!require(labdsv)){install.packages("labdsv")} ##indval
require(labdsv)

if(!require( ggplotify)){install.packages("ggplotify")}

require( ggplotify)

require(effects)

require(sjPlot)

##for figure 3
if(!require(gridExtra)){install.packages("gridExtra")}
require(gridExtra)

if(!require( ggpubr )) {install.packages("ggpubr")}
require( ggpubr )

if(!require( grid )) {install.packages("grid")}
require( grid )

if(!require( ggplot2 )) {install.packages("ggplot2")}
require( ggplot2 )

if(!require( lattice )) {install.packages("lattice")}
require( lattice )


data.genus=read_excel("/Users/naima/Projet_Diversification/Dernier_nov/native_migrates_archive/Genus&NbSample_perempo.xlsx", col_names  = TRUE)
data.genus=as.data.frame(data.genus)
row.names(data.genus)=data.genus[,1]
data.genus=data.genus[,-1]

mat_genus=data.genus[,-1]
mat_genus.1=t(mat_genus)

mat_genus.hel=decostand(mat_genus.1,"hel")

groupes=array(dim=17)
names(groupes)=rownames(mat_genus.1)

groupes[which(names(groupes) %in% c('AnimaldistalGut','Animalsurface','Animalsecretion','AnimalCorpus','Plantcorpus','Animalproximalgut','AerosolNonSaline'))]=1

groupes[which(names(groupes) %in% c('SoilNonSalin','WaterNonSalin','SedimentNonSalin','SurfaceNonSalin','Plantrhizosphere'))]=2

groupes[which(names(groupes) %in% c('SurfaceSalin','WaterSalin','Sedimentsaline','Hypersaline','Plantsurface'))]=3


iva = indval(mat_genus.hel, groupes, numitr=10000)


gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(mat_genus.hel > 0, 2, sum)[iva$pval<=0.05]


####
fidg <- data.frame(
  group = gr,
  indval = iv,
  pvalue = pv,
  freq = fr
)
fidg <- fidg[order(fidg$group, -fidg$indval), ]

fidg2=fidg
fidg2$group[fidg2$group==1]="Animal"
fidg2$group[fidg2$group==2]="Non-Saline"
fidg2$group[fidg2$group==3]="Saline"
write.xlsx(fidg2,file="indval_dfs1.xlsx")

#### corriger les noms de genera modifiés par indval
vrai_nom=read_excel("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/figure/indval_name.xlsx",col_names  = TRUE)


for(i in 1:dim(fidg2)[1])
 if (rownames(fidg2)[i] %in% vrai_nom$name_indval)  
   rownames(fidg2)[i] = vrai_nom$vrai_name[vrai_nom$name_indval==rownames(fidg2)[i]]

write.xlsx(fidg2,file="indval_dfs_vrainoms.xlsx")
#############

g1=fidg2[which(fidg2$group == "Animal"),]
g2=fidg2[which(fidg2$group =="Non-Saline"),]
g3=fidg2[which(fidg2$group == "Saline"),]

g1$genus_name=rownames(g1)
g2$genus_name=rownames(g2)
g3$genus_name=rownames(g3)

write.xlsx(g1,file="anim.xlsx")
write.xlsx(g2,file="nonsaline.xlsx")
write.xlsx(g3,file="saline.xlsx")

####heatmap with to 25 indicators for every group
indvalsummary1=rbind(g1[1:25,],g2[1:25,],g3[1:25,])

p <- ggplot(indvalsummary1, aes(group, genus_name)) + geom_tile(aes(fill = indval),
                    colour = "cadetblue1") + scale_fill_gradient(low = "cadetblue1",
                    high = "blue4")+xlab("")+ylab("")+
                    theme(axis.text.y = element_text(size=9.5),axis.text.x = element_text(size=12,face="bold"))
p
ggsave("indval_heatmap3.pdf",height = 10)

#### PCA for visualisation
rda.out=vegan::rda(mat_genus.hel)
rda.genus=rda.out
sc.rda=scores(rda.genus,scaling=1)
ssc=sc.rda$species
class(ssc)
rownames(ssc)

df.ssc=as.data.frame(ssc)

head(ssc)
dim(ssc)
write.xlsx(df.ssc,file="df.ssc.xlsx")
df.ssc$row_name=rownames(df.ssc)
head(df.ssc)
dim(df.ssc)

df.Animal=as.data.frame(df.ssc[which(rownames(ssc) %in%  g1$genus_name),])  
df.salin=as.data.frame(ssc[which(rownames(ssc) %in%  rownames(g3)),])  
df.nonsalin=as.data.frame(ssc[which(rownames(ssc) %in%  rownames(g2)),])  


write.csv(g1, "animal_genera.csv")
write.csv(g3, "salin_genera.csv")
write.csv(g2, "nonsalin_genera.csv")


########pca plot

p3 = as.ggplot( function () {
   plot(rda.genus, scaling=1, type="none",
        xlab="", ylab="",cex.lab=1.2,font.lab=1)
   title(xlab="PC1 (28.91%)",ylab="PC2 (20.52%)",line=2.5,cex.lab=1,font.lab=2)

   points(scores(rda.genus, display="sites", choices=c(1,2), scaling=1),
                    pch=2, col="black",cex=1,lwd=1)

   points(ssc[,1],ssc[,2],col="grey",cex=0.9,lwd=2)
   points(df.Animal$PC1,df.Animal$PC2,col="orange",cex=0.9,lwd=2)
   points(df.salin$PC1,df.salin$PC2,col="blue",cex=0.9,lwd=2)
   points(df.nonsalin$PC1,df.nonsalin$PC2,col="purple",cex=0.9,lwd=2)

   text(scores(rda.genus, display="sites", choices=c(1), scaling=1),
             scores(rda.genus, display="sites", choices=c(2), scaling=1),
             labels=rownames(scores(rda.genus, display="sites", scaling=1)),
             col="black", cex=0.7,font=1,offset =0.5,pos=2)
 })
save(p3,file=("pca_plot.RData"))

####plot figure 3-B


#####6 fev
### import GLMM results for animal group
load("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/glmm/animal/glmm_animal3.RData")
load("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/glmm/animal/datsc_animal.RData")
datsc1=datsc_animal
eff=allEffects(glmer.full.3)

#### import GLMM results for non saline group 
load("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/glmm/nonsaline/datsc_nonsalin.RData")
load("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/glmm/nonsaline/glmm_nonsaline3.RData")
datsc1=datsc_nonsalin
eff1=allEffects(glmer.full.status.3)
summary(glmer.full.status.3)

## import GLMM results for saline group
load("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/glmm/saline/datsc_saline.RData")
load("/Users/naima/Projet_Diversification/Dernier_nov/native_migrate7avril/glmm/saline/glmer.full.status.Sal.signif.RData")
datsc1=datsc_saline
eff2=allEffects(glmer.full.status.Sal.2)

set_theme(
  base = theme_classic(), 
  title.size=1.2,
  legend.title.face = "italic", # title font face
  legend.inside = TRUE,         # legend inside plot
  legend.color = "grey50",      # legend label color
  legend.pos = "top left",  # legend position inside plot
  axis.title.size = 1.2,
  axis.title.color = "black",
  axis.textsize = 1,
  legend.size = .7,
  legend.title.size = .8,
  geom.label.size = 3
)

plot1=plot_model(glmer.full.3,type="int",terms=c("nb_genus","status"),title="Animal",show.legend = TRUE, axis.title=c("","ASV number/genus"),colors = c("black","red"))
plot2=plot_model(glmer.full.status.3,type="int",terms=c("nb_genus","status"),title="Non saline",show.legend = FALSE,axis.title=c("Std Genus number",""),colors = c("black","red"))
plot3=plot_model(glmer.full.status.Sal.2,type="int",terms=c("nb_genus","status"),title="Saline",show.legend = FALSE,axis.title=c("",""),colors = c("black","red"))

save(plot1,file="plot1.RData")
save(plot2,file="plot2.RData")
save(plot3,file="plot3.RData")


plot4=textGrob("")
plot5=textGrob("")

lay <- rbind(c(1,1,1,5),
             c(1,1,1,6),
             c(1,1,1,5),
             c(1,1,1,5),
             c(2,3,4,5),
             c(2,3,4,5))

pdf("figure3.pdf")
grid.arrange(p3,plot1, plot2, plot3,plot4,plot5, layout_matrix = lay)
dev.off()