# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 10:22:55 2025

@author: telma
"""


import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3.csv", sep = '\t', low_memory=False)


#df_wide_filtered = df_wide.loc[df_wide['percent_alignment_score']<=80]
#df_wide_filtered['total_residue_plddt_over_80_old'] = df_wide_filtered['total_residue_plddt_over_80_old'].fillna(df_wide_filtered['total_plddt_counts_over_80_old'])
#df_wide_filtered['total_residue_plddt_over_80_new'] = df_wide_filtered['total_residue_plddt_over_80_new'].fillna(df_wide_filtered['total_plddt_counts_over_80_new'])
#df_wide_filtered['total_residue_plddt_over_80_old'] = df_wide_filtered['total_residue_plddt_over_80_old'].fillna(df_wide_filtered['residue_plddt_above_80_old'])
#df_wide_filtered['total_residue_plddt_over_80_new'] = df_wide_filtered['total_residue_plddt_over_80_new'].fillna(df_wide_filtered['residue_plddt_above_80_new'])
#df_wide_filtered.drop(columns=['total_plddt_counts_over_80_old','total_plddt_counts_over_80_new','residue_plddt_above_80_old','residue_plddt_above_80_new'], inplace=True)

#values = ['lddt', 'alntmscore', 'fident', 'protein_length',
 #         'residue_plddt_50', 'residue_plddt_60', 'residue_plddt_70', 'residue_plddt_80', 'residue_plddt_90']

#for v in values:
 #   df_wide_filtered.loc[:,'diff_'+v] = df_wide_filtered.loc[:,v+'_new'] - df_wide_filtered.loc[:,v+'_old']

#df_wide_filtered.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-FILTERED-2.tsv", sep = '\t', index=None)
df_wide_filtered['Foldseek transformed evalue new'] = np.log10(df_wide_filtered['evalue_new'])*-10
df_wide_filtered['Foldseek transformed evalue old'] = np.log10(df_wide_filtered['evalue_old'])*-10
#df_wide_filtered.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-FILTERED-2.tsv", sep = '\t', index=None)


df_wide_filtered = df_wide_filtered.rename(columns={'ptm_new':'AF3 pTM new',
                                                    'ptm_old':'AF3 pTM old',
                                                    'residue_plddt_median_new':'AF3 residue pLDDT median new',
                                                    'residue_plddt_median_old':'AF3 residue pLDDT median old',
                                                    'total_residue_plddt_over_80_new':'AF3 count of pLDDT above 80 new',
                                                    'total_residue_plddt_over_80_old':'AF3 count of pLDDT above 80 old',
                                                    'pae_mean_new':'AF3 PAE mean new',
                                                    'pae_mean_old':'AF3 PAE mean old',
                                                    'lddt_new':'Foldseek LDDT new',
                                                    'lddt_old':'Foldseek LDDT old',
                                                    'bits_new':'Foldseek bitscore new',
                                                    'bits_old':'Foldseek bitscore old',
                                                    'alntmscore_new':'Foldseek alntmscore new',
                                                    'alntmscore_old':'Foldseek alntmscore old',
                                                    'IPR_total_new':'IPR total new',
                                                    'IPR_total_old':'IPR total old',
                                                    'IPR_max_new':'IPR max new',
                                                    'IPR_max_old':'IPR max old'
                                                    })

#print(list(df_wide_filtered.columns))
#print(df_wide_filtered['above_80_new'])

AF3 = ['AF3 residue pLDDT median',
      'AF3 count of pLDDT above 80',
      'AF3 pTM',
      'AF3 PAE mean'
         ]

foldseek = ['Foldseek transformed evalue',
            'Foldseek LDDT',
            'Foldseek bitscore',
            'Foldseek alntmscore'
            ]

other = ['IPR total',
         'IPR max']

plddt_count = ['above_90',
               'above_80',
               'above_70',
               'above_60',
               'above_50']



increase_is_good = ['AF3 residue pLDDT median',
                    'AF3 pTM',
                    'AF3 count of pLDDT above 80',
                    'above 90',
                    'above 80',
                    'above 70',
                    'above 60',
                    'above 50',
                    'PSAURON in frame score',
                    'IPR total',
                    'IPR max',
                    'IPR transformed evalue min',
                    'Foldseek transformed evalue',
                    'Foldseek bitscore',
                    'Foldseek LDDT',
                    'Foldseek alntmscore'
                    ]
                    

decrease_is_good = ['AF3 PAE mean'
                    ]


def plot_scatter(list_of_names, df, figure_title, figure_height):
    # define plot size
    fig, axs = plt.subplots(nrows=len(list_of_names), ncols=3, figsize=(15, figure_height))
    plt.subplots_adjust(hspace=0.5)
    #fig.suptitle(figure_title, fontsize=20, y=0.99)
    for col_num, organism in enumerate(list(df['organism'].unique())):
        s, g, n, c = organism.split(' ')
        o = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        for row_num, ticker in enumerate(list_of_names):
            # find minimium required rows given we want 2 columns
            # define subplots
            tmp = df.loc[df['organism']==organism, [ticker+' new',ticker+' old']]
            
            if ticker in increase_is_good:
                # fill below
                axs[row_num, col_num].fill_between([0,tmp.max().max()],[0,tmp.max().max()],
                                                   facecolor="#90ee90", edgecolor="r", linewidth=0.0, alpha=0.5)
            if ticker in decrease_is_good:
                # fill above
                axs[row_num, col_num].fill_between([0,tmp.max().max()],[0,tmp.max().max()], 
                                                   np.max(tmp.max().iloc[1]),
                                                   facecolor="#90ee90", edgecolor="r", linewidth=0.0, alpha=0.5)
            
            sns.scatterplot(data=tmp[[ticker+' new',ticker+' old']],
                            x=ticker+' new', y=ticker+' old', alpha = 0.5, color='black',
                            ax = axs[row_num, col_num], legend = None)
            axs[row_num, col_num].set_xlabel(ticker + ' new', fontsize=16)
            axs[row_num, col_num].set_ylabel(ticker + ' old', fontsize=16)
            axs[row_num, col_num].annotate('n = ' + str(len(tmp[[ticker+' new',ticker+' old']].dropna())), xy=(0.7, 0.1), xycoords='axes fraction',
                        xytext=(+0.5, -0.5), textcoords='offset fontsize',
                        fontsize=16, verticalalignment='top',
                        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3.0))
        axs[0,col_num].set_title(o, fontsize=18)
    for i, ax in enumerate(axs.flat):
            ax.axline(xy1=(0, 0), slope=1, color='r', lw=2)
            alphabet = ['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)']
            ax.annotate(
                alphabet[i],
                xy=(0, 1), xycoords='axes fraction',
                xytext=(+0.5, -0.5), textcoords='offset fontsize',
                fontsize=16, fontstyle='italic', verticalalignment='top',
                bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    plt.tight_layout(pad=2)
    #plt.savefig(figure_title.replace(' ','_')+'.svg', format='svg')
    plt.savefig(figure_title.replace(' ', '_') + '.png', dpi=300)
    plt.show()
    plt.clf()


#plot_scatter(AF3, df_wide_filtered, "AlphaFold3 scores", 20)
#plot_scatter(foldseek, df_wide_filtered, "Foldseek scores", 20)
#plot_scatter(other, df_wide_filtered, "Other scores", 10)

def plot_scatter(list_of_names, df, figure_title, figure_height):
    # define plot size
    fig, axs = plt.subplots(nrows=len(list_of_names), ncols=3, figsize=(15, figure_height))
    plt.subplots_adjust(hspace=0.5)
    # fig.suptitle(figure_title, fontsize=20, y=0.99)
    for col_num, organism in enumerate(list(df['organism'].unique())):
        s, g, n, c = organism.split(' ')
        o = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        for row_num, ticker in enumerate(list_of_names):
            # find minimium required rows given we want 2 columns
            # define subplots
            tmp = df.loc[df['organism'] == organism, [ticker + '_new', ticker + '_old']].rename(columns=lambda name: name.replace('_', ' '))
            ticker = ticker.replace('_', ' ')
            if ticker in increase_is_good:
                # fill below
                axs[row_num, col_num].fill_between([0, 2000], [0, 2000],
                                                   facecolor="#90ee90", edgecolor="r", linewidth=0.0, alpha=0.5)
            if ticker in decrease_is_good:
                # fill above
                axs[row_num, col_num].fill_between([0, 2000], [0, 2000],
                                                   np.max(tmp.max().iloc[1]),
                                                   facecolor="#90ee90", edgecolor="r", linewidth=0.0, alpha=0.5)

            sns.scatterplot(data=tmp[[ticker + ' new', ticker + ' old']],
                            x=ticker + ' new', y=ticker + ' old', alpha=0.5, color='black',
                            ax=axs[row_num, col_num], legend=None)
            axs[row_num, col_num].set_xlabel(ticker + ' new', fontsize=18)
            axs[row_num, col_num].set_ylabel(ticker + ' old', fontsize=18)
            axs[row_num, col_num].annotate('n = ' + str(len(tmp[[ticker + ' new', ticker + ' old']].dropna())),
                                           xy=(0.7, 0.1), xycoords='axes fraction',
                                           xytext=(+0.5, -0.5), textcoords='offset fontsize',
                                           fontsize=16, verticalalignment='top',
                                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3.0))
        axs[0, col_num].set_title(o, fontsize=18)
    for i, ax in enumerate(axs.flat):
        ax.axline(xy1=(0, 0), slope=1, color='r', lw=2)
        alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
        ax.annotate(
            alphabet[i],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize=16, fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    plt.tight_layout(pad=2)
    #plt.savefig(figure_title.replace(' ', '_') + '.svg', format='svg')
    plt.savefig(figure_title.replace(' ', '_') + '.png', dpi=300)
    plt.show()
    plt.clf()



plot_scatter(plddt_count, df_wide_filtered, "Residue counts above thresholds", 25)
    




