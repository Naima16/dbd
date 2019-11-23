# coding=utf-8
#### produit diversity_input_test.tsv from otu_summary.emp_deblur_90bp.subset_2k.rare_5000.tsv
#### with all 2000 samples with their ASV taxonomy list
###  N.M 15 mai 2019 

import pandas as pd
import numpy  as np

def list_sample(df, index):
    return(set([x for x in df.loc[index]['list_samples'].split(',')]))

## summary file from EMP ftp site
path_otus = '~/otu_summary.emp_deblur_90bp.subset_2k.rare_5000.tsv'
df_otus = pd.read_csv(path_otus, sep='\t', index_col=0)
df_otus['samples'] = [list_sample(df_otus, index) for index in df_otus.index]

liste_seq=pd.DataFrame(columns=['sample','taxonomy'])
df_otus.samples.str
idx=df_otus.columns.get_loc("taxonomy") ##7 
idsample=df_otus.columns.get_loc("samples") ##9

for i in range(0,len(df_otus)):
 if (df_otus.iloc[i,idsample] != None):
   for element in df_otus.iloc[i,idsample] :
      indx1=liste_seq.index[liste_seq['sample']==element]       
      if ( indx1 != None) :
          liste_seq.loc[indx1,'taxonomy']=liste_seq.loc[indx1,'taxonomy'] + ","+ df_otus.iloc[i,idx]
      else :
          liste_seq=liste_seq.append({'sample':element,'taxonomy':df_otus.iloc[i,idx]},ignore_index=True)
liste_seq.to_csv("diversity_input_2k.tsv")
