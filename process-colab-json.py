# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:13:07 2024

@author: telma
"""

import json
import re
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import glob
from pathlib import Path
import pandas as pd
import scipy.stats as stats 


def progress_bar(progress, total):
    percent = (100*progress/float(total))
    bar = "█"*int(percent) + "-" * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")
    

def process_json(file, df):
    f = open(file)
    y = json.load(f)
    f.close()
    
    tp_count = len([i for i in y['plddt'] if i > 60])
    fp_count = len([i for i in y['plddt'] if i < 25])
    res_length = len(y['plddt'])
    
    d = {
    'transcript_id' : re.split("_scores",re.findall(r'TG.*', file)[0])[0],
    'protein_id' : re.split("[.-]",re.split("_scores",re.findall(r'TG.*', file)[0])[0])[0],
    'transcript_suffix' : re.split("[.-]",re.split("_scores",re.findall(r'TG.*', file)[0])[0])[1],
    'protein_length' : res_length,
    'rank' : re.split("_",re.findall(r'rank_00.', file)[0])[1],
    'ptm_model' : re.split("_",re.findall(r'TG.*', file)[0])[-3],
    'plddt' : [y['plddt']],
    'plddt_mean' : np.mean(y['plddt']),
    'plddt_std' : np.std(y['plddt']),
    'plddt_var' : np.var(y['plddt']),
    'plddt_median' : np.median(y['plddt']),
    'plddt_50' : len([v for v in y['plddt'] if (v>=50) &(v <60)]),
    'plddt_60' : len([v for v in y['plddt'] if (v >=60) & (v<70)]),
    'plddt_70' : len([v for v in y['plddt'] if (v >=70) & (v<80)]),
    'plddt_80' : len([v for v in y['plddt'] if (v >=80) & (v<90)]),
    'plddt_90' : len([v for v in y['plddt'] if v >= 90]),
    'plddt_tp_count' : tp_count,
    'plddt_fp_count' : fp_count,
    'ptm' : y['ptm'],
    'max_pae' : y['max_pae']
    }
    
    df = pd.concat([df,pd.DataFrame(d, index=[len(df)])])
    return df

def gather_data(directory):
    files = glob.glob(directory+'/*')
    df = pd.DataFrame(columns=['protein_id', 'transcript_suffix', 'rank', 'ptm_model', 'plddt_mean',
                          'plddt_std', 'plddt_var', 'plddt_median', 'ptm', 
                          'max_pae'])
    for index, file in enumerate(files):
        df = process_json(file, df)
        progress_bar(index+1,len(files))
    
    df['source'] = re.split('_',str(Path(directory).stem))[0]
    
    return df

def get_per_aa_plddt(directory):
    files = glob.glob(directory+'/*')
    plddt_per_aa = {}
    for index, file in enumerate(files):
        f = open(file)
        y = json.load(f)
        f.close()
        plddt_per_aa[Path(file).stem] = y['plddt']
        progress_bar(index+1,len(files))
    return plddt_per_aa


        

dirs = glob.glob('colabfold-results/Toxo6*_scores')

plddt_per_aa_65 = get_per_aa_plddt('colabfold-results\\Toxo65_colabfold_scores')
plddt_per_aa_66 = get_per_aa_plddt('colabfold-results\\Toxo66_colabfold_scores')


counts_65 = {'plddt_50' : len([v for v in plddt_per_aa_65 if (v>=50) &(v <60)]),
    'plddt_60' : len([v for v in plddt_per_aa_65 if (v >=60) & (v<70)]),
    'plddt_70' : len([v for v in plddt_per_aa_65 if (v >=70) & (v<80)]),
    'plddt_80' : len([v for v in plddt_per_aa_65 if (v >=80) & (v<90)]),
    'plddt_90' : len([v for v in plddt_per_aa_65 if v >= 90])}

counts_66 = {'plddt_50' : len([v for v in plddt_per_aa_66 if (v>=50) &(v <60)]),
    'plddt_60' : len([v for v in plddt_per_aa_66 if (v >=60) & (v<70)]),
    'plddt_70' : len([v for v in plddt_per_aa_66 if (v >=70) & (v<80)]),
    'plddt_80' : len([v for v in plddt_per_aa_66 if (v >=80) & (v<90)]),
    'plddt_90' : len([v for v in plddt_per_aa_66 if v >= 90])}


df_65 = gather_data('colabfold-results\\Toxo65_colabfold_scores')
df_66 = gather_data('colabfold-results\\Toxo66_colabfold_scores')
    
df_merge = df_65.merge(df_66, how='outer', suffixes = ['_65','_66'], on = ['protein_id',
 'transcript_suffix',
 'rank',
 'ptm_model',
 'plddt_mean',
 'plddt_std',
 'plddt_var',
 'plddt_median',
 'ptm',
 'max_pae',
 'transcript_id',
 'length',
 'plddt',
 'plddt_50',
 'plddt_60',
 'plddt_70',
 'plddt_80',
 'plddt_90',
 'plddt_tp_count',
 'plddt_fp_count',
 'source'])

df_merge.loc[df_merge['plddt_65'].isna(), 'plddt_65'] = df_merge.loc[df_merge['plddt_65'].isna(), 'plddt_66']
df_merge = df_merge.drop(columns='plddt_66').rename(columns={'plddt_65':'plddt'})
no_match = [i for i in df_65.protein_id if i not in list(df_66.protein_id)]
no_match += [i for i in df_66.protein_id if i not in list(df_65.protein_id)]
df_merge.loc[df_merge['protein_id'].isin(no_match), 'protein_is_in_65_and_66'] = False
df_merge.loc[df_merge['protein_is_in_65_and_66'].isna(),'protein_is_in_65_and_66'] = True

sns.set_theme(style="darkgrid")
palette={"Toxo66": '#35b778', "Toxo65": "#30678d"}
order=['Toxo65','Toxo66']
fig,ax = plt.subplots(2,2,figsize=(12, 12))
sns.scatterplot(data = df_65.merge(df_66, how='left', on=['protein_id','rank'],suffixes = ['_65','_66']), 
                x='plddt_mean_65', y = 'plddt_mean_66', hue='rank', alpha=0.8, ax=ax[0,0])
sns.scatterplot(data = df_65.merge(df_66, how='left', on=['protein_id','rank'],suffixes = ['_65','_66']), 
                x='plddt_median_65', y = 'plddt_median_66', hue='rank', alpha=0.8, ax=ax[1,0],legend=None)
sns.scatterplot(data = df_65.merge(df_66, how='left', on=['protein_id','rank'],suffixes = ['_65','_66']), 
                x='max_pae_65', y = 'max_pae_66', hue='rank', alpha=0.8, ax=ax[0,1],legend=None)
sns.scatterplot(data = df_65.merge(df_66, how='left', on=['protein_id','rank'],suffixes = ['_65','_66']), 
                x='ptm_65', y = 'ptm_66', hue='rank', alpha=0.8, ax=ax[1,1],legend=None)
plt.show()

fig,ax = plt.subplots(2,2,figsize=(10, 10))
sns.violinplot(data=df_merge[['plddt_mean','source']].melt(id_vars='source'), 
               hue='source', x='variable', y='value', split=True, inner='quart',
               palette=palette, ax=ax[0,0], legend=None, hue_order=order)
sns.violinplot(data=df_merge[['plddt_median','source']].melt(id_vars='source'), 
               hue='source', x='variable', y='value', split=True, inner='quart',
               palette=palette, ax=ax[1,0], legend=None, hue_order=order)
sns.violinplot(data=df_merge[['ptm','source']].melt(id_vars='source'), 
               hue='source', x='variable', y='value', split=True, inner='quart',
               palette=palette, ax=ax[0,1], hue_order=order)
sns.violinplot(data=df_merge[['max_pae','source']].melt(id_vars='source'), 
               hue='source', x='variable', y='value', split=True, inner='quart',
               palette=palette, ax=ax[1,1], legend=None, hue_order=order)
plt.show()


# counts 

plddt_counts = df_merge[['plddt_50','plddt_60','plddt_70', 'plddt_80', 'plddt_90','source']].melt(id_vars='source')

fig,ax = plt.subplots(1, figsize=(10, 10))
sns.barplot(data = plddt_counts.groupby(by=['source','variable'], as_index=False).sum(), 
            y='value', x='source', hue='variable', order=order, palette='magma_r')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plddt_counts_65 = df_65.merge(df_66, how='left', on=['protein_id','rank'],suffixes = ['_65','_66'])[[
    'protein_id','rank','plddt_50_65','plddt_60_65','plddt_70_65', 'plddt_80_65', 'plddt_90_65']].melt(
        id_vars=['protein_id','rank'], value_name='plddt_count_65')
plddt_counts_65['variable'] = plddt_counts_65.variable.str.replace('_65','')
plddt_counts_66 = df_65.merge(df_66, how='left', on=['protein_id','rank'],suffixes = ['_65','_66'])[[
    'protein_id','rank','plddt_50_66','plddt_60_66','plddt_70_66', 'plddt_80_66', 'plddt_90_66']].melt(
        id_vars=['protein_id','rank'], value_name='plddt_count_66')
plddt_counts_66['variable'] = plddt_counts_66.variable.str.replace('_66','')

plddt_counts = plddt_counts_65.merge(plddt_counts_66, how='outer')


grid = sns.FacetGrid(data=plddt_counts, col="variable", hue="variable", palette="magma_r",
                     col_wrap=3)
grid.map_dataframe(sns.scatterplot, x="plddt_count_65", y="plddt_count_66", marker="o",  alpha=0.7, s=3)
grid.add_legend()
grid.tight_layout()

# pointplot of counts
plddt_count_melt = plddt_counts[['protein_id', 'rank', 'variable', 'plddt_count_65', 'plddt_count_66']].melt(
    id_vars=['protein_id', 'rank', 'variable'], value_vars=['plddt_count_65', 'plddt_count_66'], 
    var_name='source', value_name='plddt_count')
sns.catplot(data=plddt_count_melt, x='source', y='plddt_count', hue='variable' , kind='point', palette="magma_r", col='rank', col_wrap=2)

# running stats

df_merge_2 = df_65.merge(df_66, on = ['protein_id', 'rank'] ,how='outer', suffixes = ['_65','_66'])[['protein_id', 'rank', 'plddt_mean_65','plddt_mean_66']].dropna()

# paired 2 test
sns.histplot(data=df_merge_2)
stats.ttest_rel(list(df_merge_2['plddt_mean_65']),list(df_merge_2['plddt_mean_66'])) 
