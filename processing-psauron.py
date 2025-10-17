# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 13:25:19 2025

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




toxo_old_psauron ="E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/TGME49-CDS-sequences/toxo65_psauron.csv"
toxo_new_psauron ="E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/TGME49-CDS-sequences/toxo66_psauron.csv"

toxo_long = "E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/non-dimer-toxo/0.8-filtered_A3_foldseek_long.tsv"
toxo_wide = "E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/non-dimer-toxo/0.8-filtered_A3_foldseek_wide.tsv"

afum_wide = "E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/AfumAf293_A3_foldseek_wide.tsv"
fgram_wide = "E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/FgramPH1-a3-and-foldseek-wide.tsv"


# aspergillus
afum_wide_df = pd.read_csv(afum_wide, sep = '\t')

afum_old_psauron = "E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Afum_old_psauron_score.csv"
afum_new_psauron ="E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Afum_new_psauron_score.csv"
afum_old_psauron = pd.read_csv(afum_old_psauron, skiprows=3)
afum_new_psauron = pd.read_csv(afum_new_psauron, skiprows=3)

afum_old_psauron['old_id'] = afum_old_psauron['description'].str.upper().str.split('-', expand =True)[0]
afum_new_psauron['new_id'] = afum_new_psauron['description'].str.upper().str.split('-', expand =True)[0]
afum_old_psauron = afum_old_psauron.rename(columns = {'in_frame_score':'in_frame_score_old',
                                                      'mean_out_of_frame_score':'mean_out_of_frame_score_old'})
afum_new_psauron = afum_new_psauron.rename(columns = {'description':'random_id','in_frame_score':'in_frame_score_new',
                                                      'mean_out_of_frame_score':'mean_out_of_frame_score_new'})

afum_wide_df = afum_wide_df.merge(afum_old_psauron[['old_id','in_frame_score_old','mean_out_of_frame_score_old']], 
                                  how='left')
afum_wide_df = afum_wide_df.merge(afum_new_psauron[['random_id','in_frame_score_new','mean_out_of_frame_score_new']], 
                                  how='left')


# fusarium
fgram_wide_df = pd.read_csv(fgram_wide, sep = '\t')

fgram_old_psauron ="E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Fgram_old_psauron_score.csv"
fgram_new_psauron ="E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/Fgram_new_psauron_score.csv"


fgram_old_psauron = pd.read_csv(fgram_old_psauron, skiprows=3)
fgram_new_psauron = pd.read_csv(fgram_new_psauron, skiprows=3)

fgram_old_psauron['old_id'] = fgram_old_psauron['description'].str.replace('T','G')
fgram_new_psauron['new_id'] = fgram_new_psauron['description'].str.replace('T','G')
fgram_old_psauron = fgram_old_psauron.rename(columns = {'in_frame_score':'in_frame_score_old',
                                                      'mean_out_of_frame_score':'mean_out_of_frame_score_old'})
fgram_new_psauron = fgram_new_psauron.rename(columns = {'in_frame_score':'in_frame_score_new',
                                                      'mean_out_of_frame_score':'mean_out_of_frame_score_new'})
fgram_new_psauron['random_id'] = fgram_new_psauron['description'].str.replace('-','_').str.upper()

fgram_wide_df = fgram_wide_df.merge(fgram_old_psauron[['old_id','in_frame_score_old','mean_out_of_frame_score_old']], 
                                  how='left')
fgram_wide_df = fgram_wide_df.merge(fgram_new_psauron[['random_id','in_frame_score_new','mean_out_of_frame_score_new']], 
                                  how='left')





# toxoplasma

toxo_old_psauron = pd.read_csv(toxo_old_psauron, skiprows=3)
toxo_new_psauron = pd.read_csv(toxo_new_psauron, skiprows=3)


toxo_old_psauron['toxo65_gene_id'] = toxo_old_psauron['description'].str.replace('-t26_1','')
toxo_old_psauron = toxo_old_psauron.rename(columns = {'in_frame_score':'in_frame_score_65',
                                                      'mean_out_of_frame_score':'mean_out_of_frame_score_65'})
toxo_new_psauron['toxo66_gene_id'] = toxo_new_psauron['description'].str.split('-', expand =True)[0]
toxo_new_psauron['toxo66_gene_id'] = toxo_new_psauron['toxo66_gene_id'].str.split('.', expand =True)[0]
toxo_new_psauron = toxo_new_psauron.rename(columns = {'in_frame_score':'in_frame_score_66',
                                                      'mean_out_of_frame_score':'mean_out_of_frame_score_66'})



toxo_wide_df = pd.read_csv(toxo_wide, sep = '\t')
toxo_wide_df = toxo_wide_df.merge(toxo_old_psauron[['toxo65_gene_id','in_frame_score_65','mean_out_of_frame_score_65']], 
                                  how='left')
toxo_wide_df = toxo_wide_df.merge(toxo_new_psauron[['toxo66_gene_id','in_frame_score_66','mean_out_of_frame_score_66']], 
                                  how='left')
toxo_wide_df['Modification_simple'] = toxo_wide_df['Modification_simple'].str.replace('modify', 'changed')

    # plot scatter plot
    mod_order=['changed', 'merge',  'split']
    
    sns.set_theme(palette='tab10', style='white')
    fig, ax = plt.subplots(3,2, figsize=(15, 15), sharex=True)
    sns.scatterplot(data=toxo_wide_df.loc[toxo_wide_df['has_clash_new']!=1,['in_frame_score_new','in_frame_score_old', 'modification_simple']], 
                    x='in_frame_score_new', y='in_frame_score_old', hue='modification_simple', alpha = 0.5,
                    ax = ax[0,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=toxo_wide_df.loc[toxo_wide_df['has_clash_new']!=1,['mean_out_of_frame_score_new','mean_out_of_frame_score_old', 'modification_simple']], 
                    x='mean_out_of_frame_score_new', y='mean_out_of_frame_score_old', hue='modification_simple', alpha = 0.5,
                    ax = ax[0,1], legend = True, hue_order=mod_order)
    sns.scatterplot(data=afum_wide_df.loc[afum_wide_df['has_clash_new']!=1,['in_frame_score_new','in_frame_score_old', 'modification']], 
                    x='in_frame_score_new', y='in_frame_score_old', hue='modification', alpha = 0.5,
                    ax = ax[1,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=afum_wide_df.loc[afum_wide_df['has_clash_new']!=1,['mean_out_of_frame_score_new','mean_out_of_frame_score_old', 'modification']], 
                    x='mean_out_of_frame_score_new', y='mean_out_of_frame_score_old', hue='modification', alpha = 0.5,
                    ax = ax[1,1], legend = None, hue_order=mod_order)
    sns.scatterplot(data=fgram_wide_df.loc[fgram_wide_df['has_clash_new']!=1,['in_frame_score_new','in_frame_score_old', 'modification']], 
                    x='in_frame_score_new', y='in_frame_score_old', hue='modification', alpha = 0.5,
                    ax = ax[2,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=fgram_wide_df.loc[fgram_wide_df['has_clash_new']!=1,['mean_out_of_frame_score_new','mean_out_of_frame_score_old', 'modification']], 
                    x='mean_out_of_frame_score_new', y='mean_out_of_frame_score_old', hue='modification', alpha = 0.5,
                    ax = ax[2,1], legend = None, hue_order=mod_order)
    sns.move_legend(ax[0,1], "upper left", bbox_to_anchor=(1, 1))
    ax[0,0].set_title('Toxoplasma gondii str. ME49', fontdict={'size': 20, 'weight': 'bold'}, x=0.3, y = 1.05)
    ax[1,0].set_title('Aspergillus fumigatus str. Af293', fontdict={'size': 20, 'weight': 'bold'}, x=0.3, y = 1.05)
    ax[2,0].set_title('Fusarium graminearum str. PH-1', fontdict={'size': 20, 'weight': 'bold'}, x=0.3, y = 1.05)
    ax[0,1].set_ylim(-0.05,1)
    ax[1,1].set_ylim(-0.05,1)
    ax[2,1].set_ylim(-0.05,1)
    for ax in ax.flat:
        ax.axline(xy1=(0, 0), slope=1, color='r', lw=2)    
    plt.tight_layout()
    plt.savefig("psauron_scatterplots.svg", format='svg')
    plt.show()
    plt.clf()
    
    # plot scatter plot
    mod_order=['changed', 'merge',  'split', 'new', 'deletion','iso_gain']
    
    sns.set_theme(palette='tab10', style='white')
    fig, ax = plt.subplots(3,2, figsize=(15, 15))
    sns.scatterplot(data=toxo_wide_df.loc[toxo_wide_df['has_clash_new']!=1,['in_frame_score_new','ph-LDDT_new', 'modification_simple']], 
                    x='in_frame_score_new', y='ph-LDDT_new', hue='modification_simple', alpha = 0.5,
                    ax = ax[0,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=toxo_wide_df.loc[toxo_wide_df['has_clash_new']!=1,['in_frame_score_old','ph-LDDT_old', 'modification_simple']], 
                    x='in_frame_score_old', y='ph-LDDT_old', hue='modification_simple', alpha = 0.5,
                    ax = ax[0,1], legend = True, hue_order=mod_order)
    sns.scatterplot(data=afum_wide_df.loc[afum_wide_df['has_clash_new']!=1,['in_frame_score_new','ph-LDDT_new', 'modification']], 
                    x='in_frame_score_new', y='ph-LDDT_new', hue='modification', alpha = 0.5,
                    ax = ax[1,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=afum_wide_df.loc[afum_wide_df['has_clash_new']!=1,['in_frame_score_old','ph-LDDT_old', 'modification']], 
                    x='in_frame_score_old', y='ph-LDDT_old', hue='modification', alpha = 0.5,
                    ax = ax[1,1], legend = None, hue_order=mod_order)
    sns.scatterplot(data=fgram_wide_df.loc[fgram_wide_df['has_clash_new']!=1,['in_frame_score_new','ph-LDDT_new', 'modification']], 
                    x='in_frame_score_new', y='ph-LDDT_new', hue='modification', alpha = 0.5,
                    ax = ax[2,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=fgram_wide_df.loc[fgram_wide_df['has_clash_new']!=1,['in_frame_score_old','ph-LDDT_old', 'modification']], 
                    x='in_frame_score_old', y='ph-LDDT_old', hue='modification', alpha = 0.5,
                    ax = ax[2,1], legend = None, hue_order=mod_order)
    sns.move_legend(ax[0,1], "upper left", bbox_to_anchor=(1, 1))
    ax[0,0].set_title('Toxoplasma gondii str. ME49', fontdict={'size': 20, 'weight': 'bold'}, x=0.3, y = 1.05)
    ax[1,0].set_title('Aspergillus fumigatus str. Af293', fontdict={'size': 20, 'weight': 'bold'}, x=0.3, y = 1.05)
    ax[2,0].set_title('Fusarium graminearum str. PH-1', fontdict={'size': 20, 'weight': 'bold'}, x=0.3, y = 1.05)
    for ax in ax.flat:
        ax.axline(xy1=(0, 0), slope=1, color='r', lw=2)    
    plt.tight_layout()
    plt.savefig("psauron_ph-LDDT_scatterplots.svg", format='svg')
    plt.show()
    plt.clf()
    
    