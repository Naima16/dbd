# coding=utf-8
import pandas as pd
import numpy  as np


def list_sample(df, index):
    #return(set([x for x in df.loc[index]['list_samples'].split(',') if x.split('.')[0] in set_study]))
    return(set([x for x in str(df.loc[index]['list_samples']).split(',')]))
  
path_otus = '~/otu_summary.emp_deblur_90bp.subset_2k.rare_5000.tsv'

df_otus = pd.read_csv(path_otus, sep='\t', index_col=0,low_memory=False)

df_otus['samples'] = [list_sample(df_otus, index) for index in df_otus.index]

col_tokeep=['sequence','samples','num_samples'] 

df_otus[col_tokeep].to_csv("seq_sample_test.tsv")

liste_seq=pd.DataFrame(columns=['sample','sequence'])

idx = 0 ##sequence

for i in range(0,len(df_otus)):
 if (df_otus.iloc[i,9] != None):
   for element in df_otus.iloc[i,9] :
      indx1=liste_seq.index[liste_seq['sample']==element]
      element_exist=0
      for k in range(len(liste_seq)) :
        
        if(liste_seq.iloc[k,0]==element) :
            indx1=liste_seq.index[liste_seq['sample']==element]
            element_exist=1
            break
            
      if ( element_exist==1) :
          liste_seq.loc[indx1,'sequence']=liste_seq.loc[indx1,'sequence']+","+df_otus.iloc[i,idx]
      else :
          liste_seq=liste_seq.append({'sample':element,'sequence':df_otus.iloc[i,idx]},ignore_index=True)

#liste_seq.to_csv("diversity_input_test.tsv")

####construction des fastas par sample
for i in range(len(liste_seq)):
   file_name=str(liste_seq.iloc[i,0])+".fa"
   ofile = open(file_name, "w")
   liste_des_seq=liste_seq.iloc[i,1].split(",") 
   for j in range(len(liste_des_seq)) :
          ofile.write(">" + liste_des_seq[j]+"\n")
          ofile.write(liste_des_seq[j]+"\n")

   ofile.close()
