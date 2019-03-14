#### construct the final array used in GLMM analyses (example OrderFamily_tab.tsv for family:order)
#### with 3 columns : sample, non focal lineage number, focal lineage sub-lineage number
#### 
####
#### NM 20 juin 2018
##########################################
##########################################

from __future__ import division
import pandas as pd
import numpy  as np


def list_otu_studies(df, index):
    return(set([x.split('.')[0] for x in df.loc[index]['list_samples'].split(',')]))

def list_taxonomy(df, index):
    return(set([x.split('__')[1] for x in df.loc[index]['taxonomy'].split(';')]))
  
## contient toutes les lineage de niveau i et lineage de niveau i-1 correspondantes    
f= open("listeOrderFamily.tsv","a")

####dbd_analys_input.py output
path_otus ='~/diversity_input_2k.tsv'
liste_seq=pd.read_csv(path_otus,sep=',')

##### compter les order_type versus classes, sur la df liste_seq à deux colonne study : entier
###### et taxonomy : chaine de caratere concatenation des taxonmy separée par "/t"
tableau_final=pd.DataFrame(columns=['sample','nb_order','nb_family'])
total_tab=pd.DataFrame(columns=['sample','order_type','nb_order','nb_family'])


for i in range(0,len(liste_seq)):
  taxonomy1=liste_seq.iloc[i,2]
  liste_taxonomy=taxonomy1.split(",") 
  liste_order_type_family=pd.DataFrame(columns=['order_type','family','sublevel_nb'])
  k=0
  break_var=0
  for j in range(0,len(liste_taxonomy)) :
     try:
       myorder_type=liste_taxonomy[j].split(";")[3].split('__')[1] ##5:genre,4:famille,3:ordre,2:class,1:phylum
     except IndexError:
       myorder_type="Unknown"
     try:
       myfamily=liste_taxonomy[j].split(";")[4].split('__')[1]  ##6:family,5:genre,4:famille,3:ordre,2:class
     except IndexError:
       myfamily="Unknown"
     ### si la taxonomy arrete a order_type on ne compte pas la family
     if (myorder_type !="Unknown" and myorder_type !="") :
       
      if (myfamily == "") :
           myfamily = "Unknown"
      if (not any (liste_order_type_family.order_type == myorder_type)): 
         liste_order_type_family=liste_order_type_family.append({'order_type':myorder_type,'family':myfamily,'sublevel_nb':1},ignore_index=True)
         
      else : 
       indx1=liste_order_type_family.index[liste_order_type_family['order_type']==myorder_type]
       liste_str=liste_order_type_family.loc[indx1,'family'].str.split(',')
       allerplusloin=1
       for kk in range(0,len(liste_str)) :
          if (myfamily in str(liste_str.iloc[kk])) :
             allerplusloin=0
             break
       if (allerplusloin==1) :
          liste_order_type_family.loc[indx1,'family']=liste_order_type_family.loc[indx1,'family'].astype(str)+","+myfamily
          ##ajouter 1 à nb_family
          liste_order_type_family.loc[indx1,'sublevel_nb'] += 1
     
  
  for k in range(0,len(liste_order_type_family)) :
       total_tab=total_tab.append({'sample':liste_seq.iloc[i,1],'order_type':liste_order_type_family.iloc[k,0],'nb_order': len(liste_order_type_family)-1,'nb_family':liste_order_type_family.iloc[k,2] },ignore_index=True)
  f.write( liste_seq.iloc[i,1] + '\n' )
  liste_order_type_family.to_csv(f)
  f.write('\n' )
  
### table containing sample, order_name, nb_order and nb_family. this table is used to fit the GLMM
total_tab.to_csv("OrderFamily_tab.tsv")
f.close()   

     
