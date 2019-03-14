### script to plot Figure 4 
### plot dbd slope in function of diversity (mean lineage number)
### NM Janv 2019

require(readxl)
require(ggplot2)
require(ggpubr)
require(effects)

### read phylum:class data [sample,phylum_name,nb_phylum,nb_class] (output of glmm_analys_input.py) 
tab_phylatype=read.table("~/phyla_tab_2k.tsv",header=T,sep=",")
tab_phylatype=as.data.frame(tab_phylatype)
tab_phylatype$phylum_type=as.factor(tab_phylatype$phylum_type)

### EMPO3 designation for samples
tab_emp=read.csv("~/sample_emplevel.tsv",header=T)
tab_emp=tab_emp[,-1]
colnames(tab_emp)=c("sample","empo_0", "empo_1" ,    "empo_2"  ,   "empo_3")

tab_emp=as.data.frame(tab_emp)
tab_emp$empo_1=as.factor(tab_emp$empo_1)
tab_emp$empo_2=as.factor(tab_emp$empo_2)
tab_emp$empo_3=as.factor(tab_emp$empo_3)

tab_phylatype_emp=merge(tab_phylatype[,-1],tab_emp[,c("sample","empo_3","empo_2","empo_1")],by="sample")
ff=tab_phylatype_emp[,c("sample","nb_phylum","empo_1","empo_3")]

####plot dbd slope in function of diversity (figure 4)

##### coefficients of the GLMM with diversity*biome as a fixed effect
##### see supllements file 1 section 4.1
data.1=read_excel("~/estimation_biomefixe.xlsx", col_names  = TRUE)

data.1$Phylum_Class = as.numeric(data.1$Phylum_Class)
data.1$SE=as.numeric(data.1$SE)
data.1$EMPO1=as.factor(data.1$EMPO1)
data.1$exp=exp(data.1$Phylum_Class)
divers=unique(ff[,c("sample","nb_phylum","empo_3","empo_1")])

r1<-with(divers, tapply(nb_phylum, empo_3, mean))

r1=as.data.frame(r1)
colnames(r1)
r1$Biome=rownames(r1)
colnames(r1)=c("diversity","Biome")

tab.total1=merge(data.1,r1,by="Biome")
tab.total1$Phylum_Class=as.numeric(tab.total1$Phylum_Class)
tab.total1$diversity=as.numeric(tab.total1$diversity)

lm1=lm(tab.total1$Phylum_Class~tab.total1$diversity)
summary(lm1)

pval=anova(lm1)$`Pr(>F)`[1]
plab1 <- sprintf('p-value: %g', pval)
radj=summary(lm1)$adj.r.squared
pr=sprintf('Adj. R2: %5.3f', radj)

p1=ggplot(tab.total1, aes(diversity, Phylum_Class,label=Biome,color=EMPO1))+
  geom_point(size = 2)+
  scale_color_manual(values=c("red", "black"))+
  geom_smooth(aes(group=1),method="lm",color="black")+
  geom_text(aes(vjust = -0.9,hjust=0.4),size=3)+
  geom_text(data = tab.total1, aes(x = +Inf, y = +Inf,label = plab1), size = 4, hjust = 1,vjust =2.5,color="black")+
  geom_text(data = tab.total1, aes(x = +Inf, y = +Inf,label = pr), size = 4, hjust = 1,vjust = 1,color="black")+
  labs(y="DBD slope",title="",x="Phylum diversity")+
  theme(legend.position = c(0.1,0.2))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(colour="black", fac="bold",size=10))
p1 
ggsave("phylum_dbd_diversity.pdf")

#pdf("DBD_diversity_phyla.pdf")
ggarrange(p1, p2,  common.legend = TRUE)
#dev.off()
ggsave("DBD_diversity_phyla.pdf")


##### genome size plot
#### GLMM output, supplementary, file1, section 5
load("~/glmersize34.RData")
#### data used in the GLMM above
load("~/datsc1.RData")

f1=allEffects(glmer.size.34)
f1=as.data.frame(f1)

tab=f1$`nb_genus:size`
tabminus2=tab[tab$size==-2,]
lm1=lm(tabminus2$fit~tabminus2$nb_genus)
summary(lm1)

tabminus1=tab[tab$size==1,]
lmplus1=lm(tabminus1$fit~tabminus1$nb_genus)
summary(lmplus1)

tabminus4=tab[tab$size==4,]
lmplus4=lm(tabminus4$fit~tabminus4$nb_genus)
summary(lmplus4)

tabminus3=tab[tab$size==3,]
lmplus3=lm(tabminus3$fit~tabminus3$nb_genus)
summary(lmplus3)

tabminus03=tab[tab$size==-0.3,]
lmsize03=lm(tabminus03$fit~tabminus03$nb_genus)
summary(lmsize03)

##### build the matrix of genome size + coef, SE (from the different lms above)
##### lms models fit (ASV_number or diversification) ~ nb_genus (diversity)
##### genome_size is z transformed

genome_size=c(-2,-0.3,1,3,4) 
est=c(0.123705,0.16821,0.20364,0.26076,0.29062)
se=c( 0.004119,0.00744,0.01070 ,0.01701,0.02079)
mat_size=as.data.frame(cbind(genome_size,est,se))

#####  z-transformed genome_size into real genome_size (atualsize)
##### load the data, before z-transformation, 
##### used in GLMM in supplementary, file 1, section 5

load("~/data_gensize.RData")
moyenne=mean(data.gensize$size)
stdV=sd(data.gensize$size)
mat_size$atualsize= (mat_size$a*stdV)+moyenne

size1=ggplot(mat_size, aes(atualsize,est))+
  geom_point(size = 2)+
  theme_bw() +
  xlab("Genome size (Mbp)")+ ylab("DBD slope")+
  geom_errorbar(width=.5, size=1,aes(ymin=est-1.96*se, ymax=est+1.96*se))
size2=size1+theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  scale_x_continuous(limits=c(0, 15))+
  theme(axis.title = element_text(face="bold",size=17))
size3=size2+geom_line(size=0.5)

#ggsave("fig4.pdf",width=13,heigh=7)

res.cor=cor.test(mat_size$atualsize,mat_size$est)

lay <- rbind(c(1,1,2))
pdf("fig4.pdf",width=13,heigh=7)
grid.arrange(p1, size3, layout_matrix = lay)
dev.off()

