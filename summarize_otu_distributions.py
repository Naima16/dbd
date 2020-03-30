import pandas as pd
import numpy as np
import biom

def make_otu_summary(path_table, path_otu_summary):
    path_otus =path_table
    table = biom.load_table(path_table)
    num_samples = len(table.ids(axis='sample'))
    # Get arrays of sample IDs and OTUs (sequences), dicts per OTU of total observations, 
    # number of samples, list of samples, and taxonomy
    otu_total_obs = {}
    otu_num_samples = {}
    otu_list_samples = {}
    samples = table.ids(axis='sample')
    otus = table.ids(axis='observation')
    print otus[1:10]
    print samples[1:10]
    
    
    for idx, cdat in enumerate(table.iter_data(axis='observation')):
        otu_total_obs[otus[idx]] = np.sum(cdat)
        otu_num_samples[otus[idx]] = np.sum(cdat > 0)
        otu_list_samples[otus[idx]] = samples[np.where(cdat > 0)[0]]
   
    otu_tax = {i: '; '.join(md['taxonomy']) for v, i, md in table.iter(axis='observation')}
    
    #print otu_tax.iloc[1:10]
    # Create Pandas DataFrame of index, sequence, total_obs, num_samples, list_samples
    df_otus = pd.DataFrame(index=np.arange(len(otus)))
    df_otus['sequence'] = [otus[i] for i in df_otus.index]
    df_otus['total_obs'] = [otu_total_obs[seq] for seq in df_otus.sequence]
    df_otus['num_samples'] = [otu_num_samples[seq] for seq in df_otus.sequence]
    df_otus['list_samples'] = [','.join(otu_list_samples[seq]) for seq in df_otus.sequence]
    df_otus['taxonomy'] = [otu_tax[seq] for seq in df_otus.sequence]
    
    
    # Add columns for total_obs_rank and num_samples_rank
    # sort by total_obs, reset index, rename index to total_obs
    df_otus = df_otus.sort_values('total_obs', ascending=False).reset_index(drop=True)
    df_otus.index.rename('total_obs_rank', inplace=True)
    # sort by num_samples, reset index, rename index to total_obs
    df_otus = df_otus.sort_values('num_samples', ascending=False).reset_index(drop=False)
    df_otus.index.rename('num_samples_rank', inplace=True)
    # keep sorted by num_samples, reset index
    df_otus = df_otus.reset_index(drop=False)
    # Add columns for total_obs_percent and num_samples_percent
    df_otus['total_obs_frac'] = df_otus['total_obs'] / df_otus['total_obs'].sum()
    df_otus['num_samples_frac'] = df_otus['num_samples'] / num_samples
    # Add 1 to the rank so they are true rank and not python-style
    df_otus['num_samples_rank'] = df_otus['num_samples_rank'] + 1
    df_otus['total_obs_rank'] = df_otus['total_obs_rank'] + 1
    # Reorder columns
    df_otus = df_otus[['sequence', 
                       'num_samples',
                       'num_samples_frac',
                       'num_samples_rank',
                       'total_obs', 
                       'total_obs_rank', 
                       'total_obs_frac',
                       'taxonomy',
                       'list_samples']]
    # Write to tsv
    df_otus.to_csv(path_otu_summary, sep='\t')

##utilisation
#path_table = '/Users/naima/deblur_raref_filt/exported/table-with-taxonomy_2.biom'
path_table = '~/data_withtaxonomy.biom'
path_otu_summary = 'otu_summary.tsv'

make_otu_summary(path_table, path_otu_summary)
