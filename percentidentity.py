##### Clustering the ASVs by cpt % identity 
##### the array produced [sample, Clusters number] is used in  
##### the percent nucleotide identity-based analysis
##### April 2018 NM

import os
import pandas as pd
import statistics
import numpy as np

## The directory containing the 2000 fasta files 
source = '/Users/naima/Diversification_rate_18mars/fasta_file_30mars' 
rapport=pd.DataFrame(columns=['percent','max','min','mean','median'])
tableau_final=pd.DataFrame(columns=['sample','nb_cluster'])
cpt = 0.6  ### 60% identity 

while cpt <= 1 :
  for root, dirs, filenames in os.walk(source):
  
    for f in filenames:
     if os.path.splitext(f)[1] == '.fa':
        fullpath = os.path.join(source, f)
        log = open(fullpath, 'r')
        cmd = 'usearch -id ' +  str(cpt) + ' -cluster_fast ' + fullpath +  ' -uc results.uc'
        print(cmd)
        os.system(cmd)
        with open('results.uc', 'r') as f1:
            lines = f1.read().splitlines()
            last_line = lines[-1]
            nb_clust = int(last_line.split()[1])+1
            sample_name=os.path.splitext(f)[0]
                
        tableau_final=tableau_final.append({'sample':sample_name,'nb_cluster':nb_clust},ignore_index=True)
        
  tableau_final.to_csv("tableau_identity_"+str(cpt)+".tsv")

  rapport=rapport.append({'percent':cpt,'max':max(tableau_final.nb_cluster),'min': min(tableau_final.nb_cluster),'mean':float("{0:.2f}".format(np.mean(tableau_final.nb_cluster.astype(np.float64)))),'median':np.median(tableau_final.nb_cluster.astype(np.float64))},ignore_index=True)
  cpt = float("{0:.2f}".format(cpt + 0.01))

rapport.to_csv("rapport_percentIdentity_range_mean_median_2.tsv")

