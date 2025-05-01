# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:25:36 2024

@author: telma
"""

import re
import numpy as np
import pandas as pd
import glob
from pathlib import Path
import matplotlib.pyplot as plt 
import seaborn as sns
import scipy.stats as stats 



def progress_bar(progress, total):
    percent = (100*progress/float(total))
    bar = "█"*int(percent) + "-" * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")


##### FIRST TEST RUN

def process_m8(directory):
    files = glob.glob(directory+'/*')
    format_output = ["query","target","alntmscore","lddt","fident","alnlen",
                    "mismatch","gapopen","qstart","qend","tstart","tend","evalue",
                    "bits","taxid","taxname","taxlineage", 'transcript_id', 'protein_id', 
                    'transcript_suffix','rank','ptm_model']
    fsdf = pd.DataFrame(columns=format_output)
    print("Processing " + str(directory))
    for index, file in enumerate(files):
        tmpdf = pd.read_csv(file, sep='\t', header=None, names = format_output)
        tmpdf = tmpdf.loc[tmpdf['lddt'] == tmpdf['lddt'].max()]
        tmpdf['source'] = str(Path(directory).stem)
        tmpdf['transcript_id'] = re.split("_r",re.findall(r'TG.*', file)[0])[0]
        tmpdf['protein_id'] = re.split("[.-]",re.split("_scores",re.findall(r'TG.*', file)[0])[0])[0]
        tmpdf['transcript_suffix'] = re.split('_r',re.split("[.-]",re.split("_scores",re.findall(r'TG.*', file)[0])[0])[1])[0]
        tmpdf['rank'] = re.split("_",re.findall(r'rank_00.', file)[0])[1]
        tmpdf['ptm_model'] = re.split("_",re.findall(r'TG.*', file)[0])[-3]
        fsdf = fsdf.merge(tmpdf, how='outer')
        progress_bar(index+1,len(files))
    print("")
    print("")
    return fsdf



dirs = glob.glob('foldseek-results/*[pdb|swissprot]')   

Toxo65_foldseek_pdb = process_m8('foldseek-results\\Toxo65_foldseek_pdb')
Toxo65_foldseek_swissprot = process_m8('foldseek-results\\Toxo65_foldseek_swissprot')
Toxo66_foldseek_pdb = process_m8('foldseek-results\\Toxo66_foldseek_pdb')
Toxo66_foldseek_swissprot = process_m8('foldseek-results\\Toxo66_foldseek_swissprot')
foldseek_df = pd.concat([Toxo65_foldseek_pdb,Toxo66_foldseek_pdb,Toxo65_foldseek_swissprot,Toxo66_foldseek_swissprot])


#### SECOND RUN WITH ALL MODELS

def process_m8(directory):
    files = glob.glob(directory)
    format_output = ["query","target","alntmscore","lddt","fident","alnlen",
                    "mismatch","gapopen","qstart","qend","tstart","tend","evalue",
                    "bits","taxid","taxname","taxlineage", 'transcript_id', 'protein_id', 
                    'transcript_suffix','ptm_model']
    fsdf = pd.DataFrame(columns=format_output)
    print("Processing " + str(directory))
    for index, file in enumerate(files):
        tmpdf = pd.read_csv(file, sep='\t', header=None, names = format_output)
        tmpdf = tmpdf.loc[tmpdf['lddt'] == tmpdf['lddt'].max()]
        #tmpdf['source'] = str(Path(glob.glob(directory)[0])).split('\\')[1]
        tmpdf['transcript_id'] = re.split("_r",re.findall(r'tg.*', file)[0])[0]
        tmpdf['transcript_id'] = tmpdf.transcript_id.str.replace('_model.cif.m8','').str.upper()
        tmpdf['protein_id'] = re.split("[.-]",re.split("_scores",re.findall(r'tg.*', file)[0])[0])[0]
        tmpdf['transcript_suffix'] = re.split('_r',re.split("[.-]",re.split("_scores",re.findall(r'tg.*', file)[0])[0])[1])[0]
        #tmpdf['rank'] = re.split("_",re.findall(r'rank_00.', file)[0])[1]
        tmpdf['ptm_model'] = re.split("_",re.findall(r'tg.*', file)[0])[-3]
        fsdf = fsdf.merge(tmpdf, how='outer')
        progress_bar(index+1,len(files))
    print("")
    print("")
    return fsdf

directory = 'foldseek-merged-results' 

Toxo65_foldseek_pdb = process_m8(directory+'\\*65_*_output-pdb\\*tg*')
Toxo65_foldseek_pdb['source'] = 'Toxo65_foldseek_pdb'

Toxo65_foldseek_swissprot = process_m8(directory+'\\*65_*_output-swissprot\\*tg*')
Toxo65_foldseek_swissprot['source'] = 'Toxo65_foldseek_swissprot'

Toxo66_foldseek_pdb = process_m8(directory+'\\*66_*_output-pdb\\*tg*')
Toxo66_foldseek_pdb['source'] = 'Toxo66_foldseek_pdb'

Toxo66_foldseek_swissprot = process_m8(directory+'\\*66_*_output-swissprot\\*tg*')
Toxo66_foldseek_swissprot['source'] = 'Toxo66_foldseek_swissprot'

foldseek_df = pd.concat([Toxo65_foldseek_pdb,Toxo66_foldseek_pdb,Toxo65_foldseek_swissprot,Toxo66_foldseek_swissprot])


foldseek_df['db_source'] = foldseek_df['source'].map({'Toxo66_foldseek_pdb': 'Toxo66', 
                           'Toxo65_foldseek_pdb' : 'Toxo65', 
                           'Toxo65_foldseek_swissprot': 'Toxo65',
                           'Toxo66_foldseek_swissprot': 'Toxo66'})

foldseek_df['protein_id'] = foldseek_df['protein_id'].str.upper()
foldseek_df['transcript_suffix'] = foldseek_df.transcript_suffix.str.replace('_model', '')



#### plotting

sns.set_theme(style="darkgrid")

sns.boxplot(data=foldseek_df, 
            x = 'source', y = 'lddt', legend=False,
            palette={"Toxo66_foldseek_pdb": '#35b778', "Toxo65_foldseek_pdb": "#30678d", 
                     "Toxo66_foldseek_swissprot": '#35b778', "Toxo65_foldseek_swissprot": "#30678d"})
plt.ylabel('max lddt')
plt.xticks(rotation=90)
plt.show()

sns.boxplot(data=foldseek_df, 
            x = 'source', y = 'alntmscore', legend=False,
            palette={"Toxo66_foldseek_pdb": '#35b778', "Toxo65_foldseek_pdb": "#30678d", 
                     "Toxo66_foldseek_swissprot": '#35b778', "Toxo65_foldseek_swissprot": "#30678d"})
plt.ylabel('alntmscore')
plt.xticks(rotation=90)
plt.show()


#df.loc[df['lddt']>0.6, 'confusion_matrix'] = 'TP'
#df.loc[df['lddt']<0.25, 'confusion_matrix'] = 'FP'
#df.loc['precision'] = 

