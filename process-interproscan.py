# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 16:49:16 2025

@author: telma
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


ipr_toxo = "E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/TGME49-interproscan/toxo-aa-INTERPROSCAN-2.tsv"
ipr_new = ["E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Aspergillus_fumigatus_Af293_2025/Afumigatus_07_05_24_INTERPROSCAN_.tsv",
          "E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Fusarium_graminearum_PH-1_2025/Fusarium_graminearum_26_07_2024-INTERPROSCAN_.tsv"]
ipr_old =["E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Aspergillus_fumigatus_Af293_2025/AfumigatusAf293_FungiDB-68_changed_proteins_INTERPROSCAN_.tsv",
            "E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Fusarium_graminearum_PH-1_2025/FgraminearumPH1_FungiDB-68_changed_proteins_INTERPROSCAN_.tsv"]

df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2.tsv", sep = '\t', low_memory=False)


# AFUM and FGRAM OLD

df_ipr_old = pd.DataFrame(columns=[0, 'IPR_max', 'IPR_total', 'IPR_evalue_min'])

for file in ipr_old:
    tmp = pd.read_csv(file, sep='\t', header=None)
    tmp = tmp.loc[tmp[11]!='-']
    # IPR length
    tmp[15] = tmp[7]-tmp[6]+1
       # 0 means an exact match 
    tmp.loc[tmp[8]=='-',8] = np.nan
    tmp[8] = tmp[8].astype('float')
    tmp[15] = tmp[15].astype('float')
    # IPR max and total
    tmp = tmp.groupby(0, as_index=False).agg(IPR_max=(15,'max'), IPR_total=(15,'sum'), IPR_evalue_min=(8,'min'))
    df_ipr_old = pd.concat([df_ipr_old, tmp])



# AFUM and FGRAM NEW
df_ipr_new = pd.DataFrame(columns=[0, 'IPR_max', 'IPR_total', 'IPR_evalue_min'])

for file in ipr_new:
    tmp = pd.read_csv(file, sep='\t', header=None)
    tmp = tmp.loc[tmp[11]!='-']
    # IPR length
    tmp[15] = tmp[7]-tmp[6]+1
       # 0 means an exact match 
    tmp.loc[tmp[8]=='-',8] = np.nan
    tmp[8] = tmp[8].astype('float')
    tmp[15] = tmp[15].astype('float')
    # IPR max and total
    tmp = tmp.groupby(0, as_index=False).agg(IPR_max=(15,'max'), IPR_total=(15,'sum'), IPR_evalue_min=(8,'min'))
    df_ipr_new = pd.concat([df_ipr_new, tmp])


# TOXO

df_ipr_toxo = pd.read_csv(ipr_toxo, sep='\t', header=None)
df_ipr_toxo = df_ipr_toxo.loc[df_ipr_toxo[11]!='-']
df_ipr_toxo[['gene', 'suffix']] = df_ipr_toxo[0].str.split('[-.]', expand=True)
df_ipr_toxo['gene'] = df_ipr_toxo['gene'].str.replace('Toxo65_','').str.replace('Toxo66_','')

# IPR length
df_ipr_toxo[15] = df_ipr_toxo[7]-df_ipr_toxo[6]+1
# IPR max and total
# df_ipr_toxo.groupby(0).agg(IPR_max=(15,'max'), IPR_total=(15,'sum'))
# 0 means an exact match 
df_ipr_toxo.loc[df_ipr_toxo[8]=='-',8] = np.nan
df_ipr_toxo[8] = df_ipr_toxo[8].astype('float')
ipr_toxo = df_ipr_toxo.groupby(0, as_index=False).agg(IPR_max=(15,'max'), IPR_total=(15,'sum'), IPR_evalue_min=(8,'min'))
ipr_toxo[['gene', 'suffix']] = ipr_toxo[0].str.split('[-.]', expand=True)
ipr_toxo['gene'] = ipr_toxo['gene'].str.replace('Toxo65_','').str.replace('Toxo66_','')



#df_wide_filtered = df_wide_filtered.merge(ipr_toxo.loc[ipr_toxo[0].str.contains('Toxo65')].rename(columns={'gene':'old_id'}), how='left')

df_ipr_old.rename(columns={0:'old_id'}, inplace=True)
df_ipr_new.rename(columns={0:'random_id'}, inplace=True)

ipr_toxo['new_id'] = ipr_toxo[0].str.upper().str.replace('-T-P1','').str.replace('TOXO66_','').str.replace('-T26_1','').str.replace('-P1','').str.split('.', expand=True)[0]
ipr_toxo['old_id'] = ipr_toxo[0].str.upper().str.replace('-T-P1','').str.replace('TOXO65_','').str.replace('-T26_1','').str.replace('-P1','').str.replace('01T','01G').str.split('.', expand=True)[0]

df_ipr_new = df_ipr_new.rename(columns={
    'IPR_max':'IPR_max_new', 'IPR_total':'IPR_total_new', 'IPR_evalue_min':'IPR_evalue_min_new'})


df_ipr_old = pd.concat([df_ipr_old, ipr_toxo.loc[ipr_toxo[0].str.contains('Toxo65'),
                                                 ['old_id','IPR_max', 'IPR_total', 'IPR_evalue_min']]]).rename(columns={
    'IPR_max':'IPR_max_old', 'IPR_total':'IPR_total_old', 'IPR_evalue_min':'IPR_evalue_min_old'})


df_wide_filtered['random_id'] = df_wide_filtered['random_id'].str.upper().str.replace('-','_')
df_ipr_new['random_id'] = df_ipr_new['random_id'].str.upper().str.replace('-','_')
df_ipr_old['old_id'] = df_ipr_old['old_id'].str.upper().str.replace('-T-P1','')
df_ipr_old['old_id'] = df_ipr_old['old_id'].str.upper().str.replace('-P1','').str.replace('_01T','_01G')


test = df_wide_filtered.merge(df_ipr_new, on=['random_id'], how='left')
test = test.merge(ipr_toxo.loc[ipr_toxo[0].str.contains('Toxo66'),
                                                 ['new_id','IPR_max', 'IPR_total', 'IPR_evalue_min']].rename(columns={'IPR_max':'IPR_max_new', 'IPR_total':'IPR_total_new', 'IPR_evalue_min':'IPR_evalue_min_new'}),
                  on = ['new_id'], how = 'left', suffixes=('','_y'))
for i in 'IPR_max_new', 'IPR_total_new', 'IPR_evalue_min_new':
    test[i].fillna(test[i+'_y'], inplace=True)
test.drop(columns=['IPR_max_new_y', 'IPR_total_new_y', 'IPR_evalue_min_new_y'], inplace=True)
test = test.merge(df_ipr_old, on='old_id', how='left')
print(df_ipr_old.loc[df_ipr_old['old_id'].str.contains('FGR')])


#for i in ['0_new', 'IPR_max_new', 'IPR_total_new', 'IPR_evalue_min_new']:
#    test[i] = test[i+'_x'].fillna(test[i+'_y'])
#    test.drop(columns=[i+'_x',i+'_y'], inplace=True)
    
#test = test.merge(df_ipr_old.rename(columns={'old_id_old':'old_id'}), on='old_id', how='left')
#test['IPR_transformed_evalue_min_new'] = np.log10(test['IPR_evalue_min_new'])*-10
#test['IPR_transformed_evalue_min_old'] = np.log10(test['IPR_evalue_min_old'])*-10

test.loc[test['IPR_evalue_min_new'] > 0.001, 'IPR_evalue_min_new'] = 10
test['IPR_evalue_min_new'].fillna(10, inplace=True)
test.loc[test['IPR_evalue_min_old'] > 0.001, 'IPR_evalue_min_old'] = 10
test['IPR_evalue_min_old'].fillna(10, inplace=True)
x = test.loc[test['IPR_evalue_min_new'] != 0, 'IPR_evalue_min_new'].min()
test.loc[test['IPR_evalue_min_new'] == 0, 'IPR_evalue_min_new'] = x.item()
x = test.loc[test['IPR_evalue_min_old'] != 0, 'IPR_evalue_min_old'].min()
test.loc[test['IPR_evalue_min_old'] == 0, 'IPR_evalue_min_old'] = x.item()

# recalculate negative logs
test['IPR_transformed_evalue_min_new'] = np.log10(test['IPR_evalue_min_new'])*-10
test['IPR_transformed_evalue_min_old'] = np.log10(test['IPR_evalue_min_old'])*-10
test['diff_IPR_min_evalue'] = test['IPR_evalue_min_new'] - test['IPR_evalue_min_old']
test['diff_IPR_transformed_min_evalue'] = test['IPR_transformed_evalue_min_new'] - test['IPR_transformed_evalue_min_old']

print(test)

test.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan.tsv", sep = '\t', index=None)

