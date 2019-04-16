#### permutation test 
#### percent identity analysis
#### kendall correlation
#### NM 2019

permutet <- function(y, x, nperm, t.vec)
{
  np=0
  tout.t=NA
  for(i in 1:nperm)
  {
    y.perm  <-  sample(y)
    temp.perm <- cor.test(x,y.perm,method="kendall",alternative="greater")
    t.perm <- temp.perm$statistic
    if (t.perm >= t.vec)
      np=np+1
  }
  return((np+1)/(nperm +1))
}

### read tables with clusters numbers for every percent identity from 75 to 100%
### these files are percentidentity.py output

tab100=read.table("~/tableau_identity_1.0.tsv",header=T,sep=",")
tab97=read.table("~/tableau_identity_0.97.tsv",header=T,sep=",")
tab95=read.table("/~/tableau_identity_0.95.tsv",header=T,sep=",")
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

t.vec=(cor.test(tab100_97$nb_100,tab100_97$nb_97,method="kendall",alternative="greater"))$statistic  
cor.out2=permutet(tab100_97$nb_97,tab100_97$nb_100,nperm=999,t.vec) ##0.001

t.vec=(cor.test(tab95_97$nb_97,tab95_97$nb_95,method="kendall",alternative="greater"))$statistic 
cor.out2=permutet(tab95_97$nb_97,tab95_97$nb_95,nperm=999,t.vec)  #0.001

t.vec=(cor.test(tab95_9$nb_95,tab95_9$nb_90,method="kendall",alternative="greater"))$statistic  
cor.out2=permutet(tab95_9$nb_95,tab95_9$nb_90,nperm=999,t.vec)  ##0.001

t.vec=(cor.test(tab85_90$nb_90,tab85_90$nb_85,method="kendall",alternative="greater"))$statistic  ##24.9
cor.out2=permutet(tab85_90$nb_90,tab85_90$nb_85,nperm=999,t.vec) #0.001

t.vec=(cor.test(tab80_85$nb_85,tab80_85$nb_80,method="kendall",alternative="greater"))$statistic  ##23.09
cor.out2=permutet(tab80_85$nb_85,tab80_85$nb_80,nperm=999,t.vec)  ##0.001

t.vec=(cor.test(tab80_75$nb_80,tab80_75$nb_75,method="kendall",alternative="greater"))$statistic  ##12.44
cor.out2=permutet(tab80_75$nb_80,tab80_75$nb_75,nperm=999,t.vec)  #0.001
