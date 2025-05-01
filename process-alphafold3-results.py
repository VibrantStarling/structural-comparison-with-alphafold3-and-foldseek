# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:24:40 2025

@author: telma
"""
# import libraries

import pandas as pd
import numpy as np
from glob import glob
import re
import sys
import argparse
import seaborn as sns
import matplotlib.pyplot as plt 
import scipy.stats as stats 
from pathlib import Path
from scipy.stats.mstats import hmean  
import json
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def progress_bar(progress, total):
    percent = (100*progress/float(total))
    bar = "█"*int(percent) + "-" * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")


# get all the data files
def get_files(directory):
    conf_json = glob(directory+"/*[0-9]_confidences.json")[0]
    sum_conf_json = glob(directory+"/*summary_confidences.json")[0]
    rank_csv = glob(directory+"/ranking_scores.csv")[0]
    cif_file = glob(directory+"/*model.cif")[0]
    return conf_json, sum_conf_json, rank_csv, cif_file


def open_json(filename):
    f = open(Path(filename))
    y = json.load(f)
    f.close()
    return y

# get PAE and atomPLDDT out of 'confidences json'
def get_confidences_json_stats(conf_json):
    # conf_json is a big list of per atom values that need summarising into smaller, managable numbers
    y = open_json(conf_json)
    # get the tp and fp count, as described in the alphafold 2 paper

    res_length = len(y['atom_plddts'])
    # get the rest of the stats
    conf_stats = {
    'atom_length' : res_length,
    'atom_plddt_mean' : np.mean(y['atom_plddts']),
    'atom_plddt_median' : np.median(y['atom_plddts']),
    'atom_plddt_min' : np.min(y['atom_plddts']),
    'atom_plddt_max' : np.max(y['atom_plddts']),
    #'pae' : y['pae'],
    'pae_mean' : np.mean(y['pae']),
    'pae_min' : np.min(y['pae']),
    'pae_max' : np.max(y['pae'])
    }
    return conf_stats

# get all stat out of 'summary confidences'
def get_summary_confidences_stats(sum_conf_json):
    y = open_json(sum_conf_json)
    # sum_conf_json contains simple summary metrics like ptm and ranking score that need extracting
    sum_conf_stats = {k:y[k] for k in ['fraction_disordered', 'has_clash', 'iptm', 'ptm', 'ranking_score'] if k in y}
    sum_conf_stats['chain_count'] = len(y['chain_iptm'])
    
    return sum_conf_stats

# get ranking scores summary?
# only useful if we want to look beyond the best model
def get_stats_from_cif(cif_file):
    dico = MMCIF2Dict(r''+cif_file)
    #l = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', 
     #         '_atom_site.label_atom_id', '_atom_site.label_alt_id', '_atom_site.label_comp_id', 
      #        '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id', 
       #       '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x', '_atom_site.Cartn_y', 
        #      '_atom_site.Cartn_z', '_atom_site.occupancy', '_atom_site.B_iso_or_equiv', 
         #     '_atom_site.auth_seq_id', '_atom_site.auth_asym_id', '_atom_site.pdbx_PDB_model_num']
    l = ['_atom_site.auth_seq_id','_atom_site.B_iso_or_equiv']
    d = dict([(key, dico[key]) for key in l])
    df = pd.DataFrame.from_dict(d, orient='index')
    df = df.transpose()
    df = df[['_atom_site.auth_seq_id','_atom_site.B_iso_or_equiv']].astype({'_atom_site.auth_seq_id': 'int','_atom_site.B_iso_or_equiv': 'float'})
    df_residue = df.groupby('_atom_site.auth_seq_id', as_index=False ).mean()
    tp_count = len([i for i in df_residue['_atom_site.B_iso_or_equiv'] if i > 60])
    fp_count = len([i for i in df_residue['_atom_site.B_iso_or_equiv'] if i < 25])
    cif_stats = {
        'peptide_length' : len(df_residue),
        'residue_plddt' : list(df_residue['_atom_site.B_iso_or_equiv']),
        'residue_plddt_mean' : df_residue['_atom_site.B_iso_or_equiv'].mean(),
        'residue_plddt_median' : np.median(df_residue['_atom_site.B_iso_or_equiv'] ),
        'residue_plddt_min' : np.min(df_residue['_atom_site.B_iso_or_equiv'] ),
        'residue_plddt_max' : np.max(df_residue['_atom_site.B_iso_or_equiv'] ),
        'residue_plddt_50' : len([v for v in df_residue['_atom_site.B_iso_or_equiv']  if (v>=50) &(v <60)]),
        'residue_plddt_60' : len([v for v in df_residue['_atom_site.B_iso_or_equiv']  if (v >=60) & (v<70)]),
        'residue_plddt_70' : len([v for v in df_residue['_atom_site.B_iso_or_equiv']  if (v >=70) & (v<80)]),
        'residue_plddt_80' : len([v for v in df_residue['_atom_site.B_iso_or_equiv']  if (v >=80) & (v<90)]),
        'residue_plddt_90' : len([v for v in df_residue['_atom_site.B_iso_or_equiv']  if v >= 90]),
        'residue_plddt_count_tp' : tp_count,
        'residue_plddt_count_fp' : fp_count}
    return cif_stats


def get_stats(directory, source):
    conf_json, sum_conf_json, rank_csv, cif_file = get_files(directory)
    
    gene = Path(directory).stem.upper()
    conf_stats = get_confidences_json_stats(conf_json)
    sum_conf_stats = get_summary_confidences_stats(sum_conf_json)
    cif_stats = get_stats_from_cif(cif_file)
    
    d = {'gene': gene,
         'source': source}
    d.update(conf_stats)
    d.update(sum_conf_stats)
    d.update(cif_stats)
    
    return d

def iter_through_subdirectories(source_dir):
    source = source_dir
    directories = glob(source_dir + '/*')
    rows = []
    failed = []
    for index, directory in enumerate(directories):
        if len(glob(directory + '/ranking_scores.csv')) == 1:
            rows.append(get_stats(directory, source))
        else:
            failed.append(Path(directory).stem)
        progress_bar(index+1,len(directories))
    #df = pd.DataFrame(rows).fillna(np.nan)
    return rows, failed

def combine_df(source_dir_list):
    failed_combined = {}
    rows_combined = []
    for source_dir in source_dir_list:
        print()
        print("\033[32m {}\033[0;0m".format("Harvesting data from " + source_dir + " "))
        rows, failed = iter_through_subdirectories(source_dir)
        failed_combined[str(source_dir)] = failed
        rows_combined += rows
    return rows_combined, failed_combined

def calculate_F1(colab_seek_merge_df):
    # calculate precioision
    colab_seek_merge_df['precision'] = colab_seek_merge_df['plddt_tp_count']/colab_seek_merge_df['protein_length']
    colab_seek_merge_df['query_length'] = colab_seek_merge_df['qend'] - colab_seek_merge_df['qstart']
    # calculate sensitivity
    count = []
    for i, row in colab_seek_merge_df.iterrows():
        count.append(sum([n > 60 for n in colab_seek_merge_df['plddt'][i][colab_seek_merge_df['qstart'][i]:colab_seek_merge_df['qend'][i]]]))
    colab_seek_merge_df['query_tp_count'] = count
    
    colab_seek_merge_df['sensitivity'] = colab_seek_merge_df['query_tp_count']/colab_seek_merge_df['query_length']
    
    # calculate the harmonic mean of precision and sensitivity (F1)
    arr = np.array(list(zip(colab_seek_merge_df['precision'],colab_seek_merge_df['sensitivity'])))
    colab_seek_merge_df['F1'] = [hmean(i) for i in arr]
    return colab_seek_merge_df


def plot(df_long, df_wide):
    # this plotting is specific for one particular comparison but may be useful reference
    sns.set_theme(style="darkgrid")
    palette={"Toxo66": '#35b778', "Toxo65": "#30678d"}
    order=['Toxo65','Toxo66']
    col_pairs = [('pae_mean_65', 'pae_mean_66'),
                 ('fraction_disordered_65', 'fraction_disordered_66'),
                 ('ptm_65', 'ptm_66'),
                 ('ranking_score_65', 'ranking_score_66'),
                 ('residue_plddt_mean_65', 'residue_plddt_mean_66'),
                 ('residue_plddt_median_65', 'residue_plddt_median_66')
                 ]

        
    for i in col_pairs:
        g = sns.pointplot(data = df_wide[['Modification_simple',i[0], i[1]]].dropna().melt(id_vars='Modification_simple'), x = 'variable', y = 'value', 
                      hue = 'Modification_simple')
        sns.move_legend(g , loc = "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.show()
        plt.clf()
        
    
    # violin plots
    fig,ax = plt.subplots(2,3,figsize=(10, 10), squeeze=False)
    sns.violinplot(data=df_long[['residue_plddt_mean','source']].melt(id_vars='source'), 
                   hue='source', x='variable', y='value', split=True, inner='quart',
                   palette=palette, ax=ax[0,0], legend=None, hue_order=order)
    sns.violinplot(data=df_long[['pae_max','source']].melt(id_vars='source'), 
                   hue='source', x='variable', y='value', split=True, inner='quart',
                   palette=palette, ax=ax[1,1], legend=None, hue_order=order)
    sns.violinplot(data=df_long[['ptm','source']].melt(id_vars='source'), 
                   hue='source', x='variable', y='value', split=True, inner='quart',
                   palette=palette, ax=ax[0,2], hue_order=order)
    sns.violinplot(data=df_long[['residue_plddt_median','source']].melt(id_vars='source'), 
                   hue='source', x='variable', y='value', split=True, inner='quart',
                   palette=palette, ax=ax[0,1], legend=None, hue_order=order)
    sns.violinplot(data=df_long[['fraction_disordered','source']].melt(id_vars='source'), 
                   hue='source', x='variable', y='value', split=True, inner='quart',
                   palette=palette, ax=ax[1,2], legend=None, hue_order=order)
    sns.violinplot(data=df_long[['ranking_score','source']].melt(id_vars='source'), 
                   hue='source', x='variable', y='value', split=True, inner='quart',
                   palette=palette, ax=ax[1,0], legend=None, hue_order=order)
    sns.move_legend(ax[0,2], "upper left", bbox_to_anchor=(1, 1))
    plt.show()
    plt.clf()
    
    # point plots of counts of plddt thresholds
    plddt_count_melt = df_long[['gene', 'source','residue_plddt_50', 'residue_plddt_60', 'residue_plddt_70',
    'residue_plddt_80', 'residue_plddt_90']].melt(
        id_vars=['gene', 'source'], value_vars=['residue_plddt_50', 'residue_plddt_60', 'residue_plddt_70',
        'residue_plddt_80', 'residue_plddt_90'], 
        var_name='plddt_group', value_name='plddt_count')
    fig,ax = plt.subplots(figsize=(10, 10))
    sns.pointplot(data=plddt_count_melt, x='source', y='plddt_count', hue='plddt_group' , palette="magma_r")
    plt.show()
    plt.clf()
    

    # plot kde plot
    mod_order=['modify', 'merge', 'split', 'new']
    
    sns.set_theme(palette='tab10', style='white')
    fig,ax = plt.subplots(2,3, figsize=(15, 15))
    sns.scatterplot(data=df_wide[['residue_plddt_mean_66','residue_plddt_mean_65', 'Modification_simple']], 
                    x='residue_plddt_mean_66', y='residue_plddt_mean_65', hue='Modification_simple', alpha = 0.5,
                    ax = ax[0,0], legend = None, hue_order=mod_order)
    sns.scatterplot(data=df_wide[['residue_plddt_median_66','residue_plddt_median_65', 'Modification_simple']], 
                    x='residue_plddt_median_66', y='residue_plddt_median_65', hue='Modification_simple',alpha = 0.5,
                    ax = ax[0,1], legend = None, hue_order=mod_order)
    sns.scatterplot(data=df_wide[['ptm_66','ptm_65', 'Modification_simple']], 
                    x='ptm_66', y='ptm_65', hue='Modification_simple', alpha = 0.5,
                    ax = ax[0,2], hue_order=mod_order)
    sns.scatterplot(data=df_wide[['pae_mean_66','pae_mean_65', 'Modification_simple']], 
                    x='pae_mean_66', y='pae_mean_65', hue='Modification_simple', alpha = 0.5,
                    ax = ax[1,1], legend = None, hue_order=mod_order)
    sns.scatterplot(data=df_wide[['fraction_disordered_66','fraction_disordered_65', 'Modification_simple']], 
                    x='fraction_disordered_66', y='fraction_disordered_65', hue='Modification_simple', alpha = 0.5,
                    ax = ax[1,2], legend = None, hue_order=mod_order)
    sns.scatterplot(data=df_wide.loc[df_wide['has_clash_66']!=1,['ranking_score_66','ranking_score_65', 'Modification_simple']], 
                    x='ranking_score_66', y='ranking_score_65', hue='Modification_simple', alpha = 0.5,
                    ax = ax[1,0], legend = None, hue_order=mod_order)
    for ax in ax.flat:
        ax.axline(xy1=(0, 0), slope=1, color='r', lw=2)
    plt.tight_layout()
    plt.show()
    plt.clf()
    
    
    ## plot the differences between 65 and 66
    
    df = df_wide[['diff_residue_plddt_mean','diff_ptm','diff_iptm', 'diff_fraction_disordered','Modification']].dropna()
    
    name_map = {'modify, translational start shifted downstream': 'modify, translational start',
                   'modify, translational start shifted upstream' : 'modify, translational start',
                   'merge, isoform, annotation contributed by Nicolas Dos Santos Pacheco and Albert Tell I Puig':'merge, isoform',
                   'extended to start of gap, modify': 'modify, extended to start of gap',
                   'merge': 'merge', 
                   'modify': 'modify', 
                   'modify, isoform' : 'modify, isoform', 
                   'modify, different frame': 'modify, frame change',
                   'modify, incorrect frame': 'modify, frame change', 
                   'split': 'split', 
                   'merge, isoform': 'merge, isoform'}
    

    
    g = sns.catplot(data = df[['diff_residue_plddt_mean','diff_ptm','diff_iptm', 'diff_fraction_disordered','Modification']].melt(id_vars=['Modification']),
                y = 'value' , col = 'variable', hue = 'Modification', palette = 'tab10',kind = 'box', col_wrap=2, height=4, aspect=.6)
    plt.tight_layout()
    sns.move_legend(g, loc="upper left", bbox_to_anchor=(0, 0))
    plt.show()
    plt.clf()
    
    sns.catplot(data = df[['diff_residue_plddt_mean','Modification']].melt(id_vars=['Modification']), y = 'value' , col = 'variable', hue = 'Modification', kind = 'box', height=4, aspect=1)

    plt.show()
    plt.clf()
    
    return

def load_test_files():
    structural_csv = 'E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/0-structural-changes-metadata/TgondiiME49_structural_changes_aa.csv'
    df_structure = pd.read_csv(structural_csv)
    df_long = pd.read_csv('E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/all_and_merged_protein_alphafold3_summary_data_LONG.csv')
    df_wide = pd.read_csv('E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/all_and_merged_protein_alphafold3_summary_data_WIDE.csv')
    a3_seek_merge_df = pd.read_csv('E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/Alphafold3-and-foldseek-results-all-and-merged-LONG.csv')
    a3_seek_merge_df = pd.read_csv('E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/Alphafold3-and-foldseek-results-all-and-merged-with-F1-LONG.csv')
    colab_seek_merge_df = pd.read_csv('E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/1-small-scale-test/colabfold_foldseek_merged_df.csv')
    all_and_merged_foldseek = pd.read_csv('E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/2-all-proteins/FINAL_all_and_merged_foldseek_summary_results.csv')
    return

# structural_csv, df_long, df_wide = load_test_files()

def make_wide_df(df_long, structural_csv):
    df_structure = pd.read_csv(structural_csv)
    df_long['gene'] = df_long['gene'].str.replace('-T26_1','')
    df_long['protein_length'] = [len(x) for x in df_long['residue_plddt']]
    df_wide = df_structure[['Gene ID', 'Previous ID', 'Modification']].rename(columns = {'Gene ID': 'toxo66_gene_id','Previous ID' : 'toxo65_gene_id'})
    df_wide = df_wide.merge(df_long.loc[df_long['source'].str.contains('66')].rename(columns = {'protein_id': 'toxo66_gene_id'}), 
                            on = 'toxo66_gene_id', how = 'outer')
    df_wide = df_wide.merge(df_long.loc[df_long['source'].str.contains('65')].rename(columns = {'protein_id': 'toxo65_gene_id'}), how = 'outer', 
                            left_on = ['toxo65_gene_id'], 
                            right_on =['toxo65_gene_id'], suffixes = ('_66', '_65'))
    df_wide.insert(0, 'pair_id', list(zip( df_wide['toxo66_gene_id'], df_wide['toxo65_gene_id'])))
    #df_wide.merge(df_long.loc[df_long['source']=='Toxo66'].rename(columns = {'protein_id': 'toxo66_gene_id'}), on = 'toxo66_gene_id', how = 'outer')
    df_wide['diff_atom_plddt'] = df_wide.atom_plddt_mean_66 - df_wide.atom_plddt_mean_65
    df_wide['diff_residue_plddt_mean'] = df_wide.residue_plddt_mean_66 - df_wide.residue_plddt_mean_65
    df_wide['diff_ptm'] = df_wide.ptm_66 - df_wide.ptm_65
    df_wide['diff_pae_max'] = df_wide.pae_max_66 - df_wide.pae_max_65
    df_wide['diff_fraction_disordered'] = df_wide.fraction_disordered_66 - df_wide.fraction_disordered_65
    df_wide['diff_residue_plddt_50'] = df_wide.residue_plddt_50_66 - df_wide.residue_plddt_50_65
    df_wide['diff_residue_plddt_60'] = df_wide.residue_plddt_60_66 - df_wide.residue_plddt_60_65
    df_wide['diff_residue_plddt_70'] = df_wide.residue_plddt_70_66 - df_wide.residue_plddt_70_65
    df_wide['diff_residue_plddt_80'] = df_wide.residue_plddt_80_66 - df_wide.residue_plddt_80_65
    df_wide['diff_residue_plddt_90'] = df_wide.residue_plddt_90_66 - df_wide.residue_plddt_90_65
    df_wide['total_plddt_counts_over_80_65'] = df_wide['residue_plddt_80_65'] + df_wide['residue_plddt_90_65']
    df_wide['total_plddt_counts_over_80_66'] = df_wide['residue_plddt_80_66'] + df_wide['residue_plddt_90_66']
    df_wide['diff_total_plddt_counts_over_80'] = df_wide['total_plddt_counts_over_80_66'] -df_wide['total_plddt_counts_over_80_65']
    df_wide['diff_protein_length'] = df_wide['protein_length_66']  -df_wide['protein_length_65']
    df_wide['Modification'] = df_wide['Modification'].str.replace(', annotation contributed by Nicolas Dos Santos Pacheco and Albert Tell I Puig','')
    df_wide['Modification'] = df_wide['Modification'].str.replace('extended to start of gap, modify','modify, translational start shifted upstream')

    name_map = {'new, incorrect strand':'new', 'new': 'new',
        'modify, translational start shifted downstream': 'modify',
                   'modify, translational start shifted upstream' : 'modify',
                   'merge, isoform, annotation contributed by Nicolas Dos Santos Pacheco and Albert Tell I Puig':'merge',
                   'extended to start of gap, modify': 'modify',
                   'merge': 'merge', 
                   'modify': 'modify', 
                   'modify, isoform' : 'modify', 
                   'modify, different frame': 'modify',
                   'modify, incorrect frame': 'modify', 
                   'split': 'split', 
                   'merge, isoform': 'merge'}
    
    df_wide['Modification_simple'] = df_wide['Modification'].map(name_map)
    return df_wide

def plot_pca(df_wide, n, label_col):
     df= df_wide[['diff_residue_plddt_mean','diff_total_plddt_counts_over_80', 'diff_ptm','diff_fraction_disordered','Modification']].dropna()
     X =  StandardScaler().fit_transform(np.array(df[['diff_total_plddt_counts_over_80','diff_residue_plddt_mean','diff_ptm','diff_fraction_disordered']]))
     y = np.array(df['Modification']).reshape(-1, 1)
     
     
     df2 = df_long[['atom_plddt_mean','residue_plddt_mean','ptm', 'ranking_score', 'fraction_disordered', 'source']].dropna()
     X =  StandardScaler().fit_transform(np.array(df2[['atom_plddt_mean','residue_plddt_mean','ptm','iptm', 'ranking_score', 'fraction_disordered']]))
     y = np.array(df2['source']).reshape(-1, 1)
     
     name_map = {'merge': 8, 
                 'modify, translational start shifted downstream': 2,
                'modify': 1, 
                'modify, translational start shifted upstream' : 3,
                'modify, isoform' : 4, 
                'modify, different frame':5,
                'modify, incorrect frame': 6, 
                'split': 7, 
                'merge, isoform': 9,
                'merge, isoform, annotation contributed by Nicolas Dos Santos Pacheco and Albert Tell I Puig':9,
                'extended to start of gap, modify': 10}
     
     pca = PCA(n_components=3)
     pca.fit(X)
     X_r = pca.fit(X).transform(X)


     pca_df = pd.DataFrame(np.concatenate([X_r, y], axis=1))
     
     fig = plt.figure(figsize=(8, 6))
     ax = fig.add_subplot(111, projection = '3d')
     ax.scatter(pca_df[0], pca_df[1], pca_df[2],c=pca_df[3].map(name_map), cmap='Set1', marker='o')
     plt.legend(name_map)
     
     fig = plt.figure(figsize=(8, 6))
     ax = fig.add_subplot(111, projection = '3d')
     ax.scatter(pca_df[0], pca_df[1], pca_df[2],c=pca_df[3].map({"Toxo66": '#35b778', "Toxo65": "#30678d"}), marker='o')
     ax.legend(ax.legend_elements(),
                 loc="lower left", title="Source")
     
     pca_df[3] = pca_df[3].map(name_map)

     sns.pairplot(data = pca_df, hue=3, aspect=2)

def main():
    parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
Harvest alphafold3 data in a given list of directories
-------------------------------------------------------------
     ''',
     epilog="written by Helen Rebecca Davison") 
    parser.add_argument('-d', '--source_dir_list', \
                    nargs='+', action='append',
                    help='<Required> Set flag',
                    required=True)
    parser.add_argument('-p','--plot',  action='store', nargs='*',\
                    help="If true, plots will be made.")
    parser.add_argument('-o','--output',\
                        help="Name for the output files")
    args = parser.parse_args()
    
    source_dir_list = args.source_dir_list[0]
    print(source_dir_list)
    plot = args.plot
    out = args.out
    
    sys.tracebacklimit = 0
    
    for source_dir in source_dir_list:
        if Path(source_dir).exists():
            if Path(source_dir).is_dir():
                print()
            else:
                raise Exception(source_dir + " is not a diresctory!")
        else:
            raise FileNotFoundError(source_dir + " can't be found!")
    
    structural_csv = 'E:/Dropbox/1-Veupath-files/PROJECTS/Alphafold-toxo-test/0-structural-changes-metadata/TgondiiME49_structural_changes_aa.csv'
    df_structure = pd.read_csv(structural_csv)
    
    rows_combined, failed_combined = combine_df(source_dir_list)
    df_long = pd.DataFrame(rows_combined).fillna(np.nan)
    df_long.insert(1, 'protein_id', df_long['gene'].str.split('-', n=1, expand=True)[0])
    df_wide = make_wide_df(df_long, structural_csv)

    df_long.to_csv(out+'_long.csv', index=None, header=True)
    
    return

main()