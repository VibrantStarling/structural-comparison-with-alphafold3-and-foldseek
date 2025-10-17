import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

#df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-2.csv", sep = '\t', low_memory=False)
df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3.csv", sep = '\t', low_memory=False)


names = ['F1','alntmscore','alnlen','residue_plddt_mean', 'residue_plddt_median','pae_mean','ptm', 'fraction_disordered', 'mean_disorder_metapredict3']
def box_disorder_plot():
    fig,ax = plt.subplots(2,1, figsize=(16, 8),)
    tmp_a3 = df_wide_filtered[['fraction_disordered_old','fraction_disordered_new','organism']].melt(id_vars='organism')
    tmp_meta = df_wide_filtered[['mean_disorder_metapredict3_old', 'mean_disorder_metapredict3_new', 'organism']].melt(
        id_vars='organism')
    tmp_a3['variable'] = tmp_a3['variable'].str.replace('fraction_disordered_old','Fraction disordered, old').str.replace('fraction_disordered_new','Fraction disordered, new')
    tmp_a3 = tmp_a3.rename(columns={'organism':'Organism', 'variable' : 'Variable','value':'Proportion'})
    tmp_meta['variable'] = tmp_meta['variable'].str.replace('mean_disorder_metapredict3_old',
                                                        'Metapredict3, old').str.replace(
        'mean_disorder_metapredict3_new', 'Metapredict3, new')
    tmp_meta = tmp_meta.rename(columns={'organism': 'Organism', 'variable': 'Variable', 'value': 'Proportion'})
    sns.boxplot(data=tmp_a3,
                   x='Organism', y='Proportion',hue = 'Variable',
                   palette={"Fraction disordered, new": '#35b778', "Fraction disordered, old": "#30678d"}, ax=ax[0])
    sns.boxplot(data=tmp_meta,
                   x='Organism', y='Proportion', hue='Variable',
                   palette={"Metapredict3, new": '#35b778', "Metapredict3, old": "#30678d"}, ax=ax[1])
    labels = []
    for name in df_wide_filtered.organism.unique():
        s, g , n, c  = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c)
    ax[0].set_xticklabels(labels)
    ax[1].set_xticklabels(labels)
    sns.move_legend(ax[1], 'upper left', bbox_to_anchor=(1,1))
    sns.move_legend(ax[0], 'upper left', bbox_to_anchor=(1, 1))
    plt.show()

#box_disorder_plot()

def box_disorder_plot():
    sns.set(font_scale=1.2)
    sns.set_style('white')
    fig,ax = plt.subplots(3,1, figsize=(13, 13), sharex=True)

    tmp_prot = df_wide_filtered[['protein_length_old','protein_length_new','organism']].melt(id_vars='organism')
    #tmp_prot['value'] = np.log10(tmp_prot['value'])

    tmp_prot['variable'] = tmp_prot['variable'].str.replace('protein_length_old','Protein length, old').str.replace('protein_length_new','Protein length, new')
    tmp_prot = tmp_prot.rename(columns={'organism':'Organism', 'variable' : 'Variable','value':'Protein length (log)'})

    tmp_meta = df_wide_filtered[['mean_disorder_metapredict3_old', 'mean_disorder_metapredict3_new', 'organism']].melt(
        id_vars='organism')
    tmp_meta['variable'] = tmp_meta['variable'].str.replace('mean_disorder_metapredict3_old',
                                                        'Metapredict3, old').str.replace(
        'mean_disorder_metapredict3_new', 'Metapredict3, new')
    tmp_meta = tmp_meta.rename(columns={'organism': 'Organism', 'variable': 'Variable', 'value': 'Mean Metapredict3 score'})

    tmp_cds = df_wide_filtered[['CDS_counts_old', 'CDS_counts_new', 'organism']].melt(
        id_vars='organism')
    tmp_cds['variable'] = tmp_cds['variable'].str.replace('CDS_counts_old', 'CDS count, old').str.replace(
        'CDS_counts_new', 'CDS count, new')
    tmp_cds = tmp_cds.rename(
        columns={'organism': 'Organism', 'variable': 'Variable', 'value': 'Number of CDS per protein'})

    sns.boxplot(data=tmp_prot,
                   x='Organism', y='Protein length (log)',hue = 'Variable', log_scale=True,
                   palette={"Protein length, new": '#35b778', "Protein length, old": "#30678d"}, ax=ax[0], legend=False)
    sns.boxplot(data=tmp_meta,
                   x='Organism', y='Mean Metapredict3 score', hue='Variable',
                   palette={"Metapredict3, new": '#35b778', "Metapredict3, old": "#30678d"}, ax=ax[1], legend=False)
    sns.boxplot(data=tmp_cds,
                x='Organism', y='Number of CDS per protein', hue='Variable',
                palette={"CDS count, new": '#35b778', "CDS count, old": "#30678d"}, ax=ax[2], legend=False)
    labels = []
    for name in df_wide_filtered.organism.unique():
        s, g , n, c  = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c)
    alphabet = ['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)']
    for i in list(range(0,3)):
        ax[i].set_xticklabels(labels)
        #sns.move_legend(ax[i], 'upper left', bbox_to_anchor=(1, 1))
        ax[i].annotate(
            alphabet[i],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    new_patch = mpatches.Patch(color='#35b778', label='New models')
    old_patch = mpatches.Patch(color='#30678d', label='Old models')
    ax[0].legend(handles=[new_patch, old_patch], loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig('metadata-boxplots.png', dpi=300)
    plt.show()

#box_disorder_plot()

def scatter_disorder_plot():
    sns.set_palette('colorblind')
    fig,ax = plt.subplots(1,2, figsize=(16, 8))
    sns.scatterplot(data=df_wide_filtered[['mean_disorder_metapredict3_old','ptm_old','organism']].rename(columns={'mean_disorder_metapredict3_old':'Metapredict3, old',
                                                                                                                   'ptm_old':'pTM, old',
                                                                                                                   'organism':'Organism'}),
                   x='Metapredict3, old', y='pTM, old', hue = 'Organism', ax=ax[0],
                    hue_order=list(df_wide_filtered.organism.unique()), legend=None, alpha=0.5)
    sns.scatterplot(data=df_wide_filtered[['mean_disorder_metapredict3_new', 'ptm_new', 'organism']].rename(columns={'mean_disorder_metapredict3_new':'Metapredict3, new',
                                                                                                                   'ptm_new':'pTM, new',
                                                                                                                   'organism':'Organism'}),
                    x='Metapredict3, new', y='pTM, new', hue='Organism', ax=ax[1], alpha=0.5,
                    hue_order=list(df_wide_filtered.organism.unique()))
    labels=[]
    for name in df_wide_filtered.organism.unique():
        s, g , n, c  = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c)
    alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
    ax[0].annotate(
            alphabet[0],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top', fontfamily='sanserif',
            bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    ax[1].annotate(
        alphabet[1],
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', fontstyle='italic', verticalalignment='top', fontfamily='sanserif',
        bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    plt.tight_layout()
    sns.move_legend(ax[1], 'upper left', bbox_to_anchor=(1, 1), labels=labels)
    plt.show()

#scatter_disorder_plot()


# this need esmfold results in it
#
def plot_pointplot(df):
    df_melt = df[['organism','residue_plddt_50_new', 'residue_plddt_50_old',
    'residue_plddt_60_new',	'residue_plddt_60_old',
    'residue_plddt_70_new',	'residue_plddt_70_old',
    'residue_plddt_80_new',	'residue_plddt_80_old',
    'residue_plddt_90_new',	'residue_plddt_90_old',
                  'esmfold_residue_plddt_50_new', 'esmfold_residue_plddt_50_old',
                  'esmfold_residue_plddt_60_new', 'esmfold_residue_plddt_60_old',
                  'esmfold_residue_plddt_70_new', 'esmfold_residue_plddt_70_old',
                  'esmfold_residue_plddt_80_new', 'esmfold_residue_plddt_80_old',
                  'esmfold_residue_plddt_90_new', 'esmfold_residue_plddt_90_old']].melt(id_vars='organism')

    df_melt.loc[df_melt['variable'].str.contains('_new'), 'source'] = 'New'
    df_melt.loc[df_melt['variable'].str.contains('_old'), 'source'] = 'Old'
    df_melt.loc[df_melt['variable'].str.contains('esmfold_'), 'software'] = 'ESMfold'
    df_melt.loc[~df_melt['variable'].str.contains('esmfold_'), 'software'] = 'Alphafold3'
    df_melt['variable'] = df_melt['variable'].str.replace('_new','').str.replace('_old','')
    df_melt['variable'] = df_melt['variable'].str.replace('esmfold_residue_plddt_50', '≥50 to <60')
    df_melt['variable'] = df_melt['variable'].str.replace('esmfold_residue_plddt_60', '≥60 to <70')
    df_melt['variable'] = df_melt['variable'].str.replace('esmfold_residue_plddt_70', '≥70 to <80')
    df_melt['variable'] = df_melt['variable'].str.replace('esmfold_residue_plddt_80', '≥80 to <90')
    df_melt['variable'] = df_melt['variable'].str.replace('esmfold_residue_plddt_90', '≥90')
    df_melt.rename(columns={'organism':'Organism',
                            'source':'Source',
                            'variable':'Threshold',
                            'value':'Number of Residues'})
    g = sns.catplot(df_melt, x='source', y='value', hue='variable', order=['old','new'], col='software', kind='point', palette='magma')
    for label, ax in g.items():
        # Use Axes.annotate to put the label
        # - at the top left corner (axes fraction (0, 1)),
        # - offset half-a-fontsize right and half-a-fontsize down
        #   (offset fontsize (+0.5, -0.5)),
        # i.e. just inside the axes.
        ax.annotate(
            label,
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
    plt.show()
    return

#plot_pointplot(df_wide_filtered)

#print(df_wide_filtered.groupby('organism').count())

def pairplot_best(df):
    #    tmp = df[['organism','total_residue_plddt_over_80_new',
    #              'bits_new', 'IPR_total_new']].rename(columns={'organism':'Organism',
    #                                           'total_residue_plddt_over_80_new':'AF3 count of pLDDT above 80 new',
    #                                           'bits_new':'Foldseek bitscore new',
    #                                           'IPR_total_new':'IPR total new',
    #                                           })
    #    sns.pairplot(tmp, hue='Organism', corner=True)
    #    plt.savefig('pairplot-best-scores.png', dpi=300)
    #    plt.show()

    tmp = df[['organism', 'total_residue_plddt_over_80_new',
              'bits_new', 'IPR_total_new']].rename(columns={'organism': 'Organism',
                                                            'total_residue_plddt_over_80_new': 'AF3 count of pLDDT above 80 new',
                                                            'bits_new': 'Foldseek bitscore new',
                                                            'IPR_total_new': 'IPR total new',
                                                            })
    x_vars = ["Foldseek bitscore new", "IPR total new"]
    y_vars = ["AF3 count of pLDDT above 80 new"]
    g = sns.PairGrid(tmp, diag_sharey=False, corner=True, height=4, hue = 'Organism')
    g.map_lower(sns.scatterplot, alpha = 0.3, s = 5)
    g.map_diag(sns.histplot, multiple="stack", element="step")
    index=0
    for ax in g.axes.flat:
        if ax is not None:
            alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
            ax.annotate(
                alphabet[index],
                xy=(0, 1), xycoords='axes fraction',
                xytext=(+0.5, -0.5), textcoords='offset fontsize',
                fontsize='large', fontstyle='italic', verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8,edgecolor='none', pad=3.0))
            index+=1
    plt.tight_layout()
    plt.savefig('pairplot-best-scores-and-disorder-scatter.png', dpi=300)
    plt.show()
    #fig, axs = plt.subplots(1,2, figsize=(10,6))
    #sns.regplot(tmp, y = "AF3 count of pLDDT above 80 new", x = "Foldseek bitscore new",
    #            scatter_kws={'alpha':0.15, 's':5},
    #            ax = axs[0])
    #sns.regplot(tmp, y="AF3 count of pLDDT above 80 new", x="IPR total new",
    #            scatter_kws={'alpha': 0.15, 's': 5},
    #            ax=axs[1])
    #for i, x in enumerate(x_vars):
        #r, _ = stats.pearsonr(tmp[x].fillna(0), tmp["AF3 count of pLDDT above 80 new"])
        #axs[i].annotate("r = {:.2f}".format(r),
         #           xy=(0.7, 0.1), xycoords=axs[i].transAxes,
          #          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3.0))
     #   alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
      #  axs[i].annotate(
       #     alphabet[i],
        #    xy=(0, 1), xycoords='axes fraction',
         #   xytext=(+0.5, -0.5), textcoords='offset fontsize',
          #  fontsize='large', fontstyle='italic', verticalalignment='top',
           # bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    #plt.tight_layout()
    #plt.savefig('pairplot-best-scores-and-disorder-lowess.png', dpi=300)
    #plt.show()
    return


def af3_paired_disorder_plot(df):
    tmp = df[['organism', 'total_residue_plddt_over_80_new',
              'ptm_new', 'residue_plddt_median_new', 'pae_mean_new', 'mean_disorder_metapredict3_new']].rename(columns={'ptm_new':'AF3 pTM new',
                                                    'residue_plddt_median_new':'AF3 residue pLDDT median new',
                                                    'total_residue_plddt_over_80_new':'AF3 count of pLDDT above 80 new',
                                                    'pae_mean_new':'AF3 PAE mean new',
                                                    'mean_disorder_metapredict3_new':'mean metapredict3 new'})
    sns.set_context('talk')
    def corrfunc(x, y, **kws):
        r, _ = stats.pearsonr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.2f}".format(r),
                    xy=(0.7, 0.1), xycoords=ax.transAxes,bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3.0))

    g = sns.PairGrid(tmp, diag_sharey=False, corner=True, height=4)
    g.map_lower(sns.regplot, lowess=True, scatter_kws={'alpha': 0.15, 's': 5})
    g.map_lower(corrfunc, hue=None)
    g.map_diag(sns.histplot, multiple="stack", element="step")
    index=0
    for ax in g.axes.flat:
        if ax is not None:
            alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
            ax.annotate(
                alphabet[index],
                xy=(0, 1), xycoords='axes fraction',
                xytext=(+0.5, -0.5), textcoords='offset fontsize',
                fontsize='large', fontstyle='italic', verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8,edgecolor='none', pad=3.0))
            index+=1
    plt.savefig('pairplot-AF3-scores-and-disorder-lowess.png', dpi=300)
    plt.show()
    return

def pairplot_maxabs_diffs():
    df = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3-maxabs.tsv",
                     sep='\t')
    tmp = df[['AF3 total residue plddt over 80 diff.',
              'Foldseek bitscore diff.', 'IPR domain total length diff.']].rename(columns={'organism': 'Organism'})
    tmp['Organism'] = df_wide_filtered['organism']
    labels = []
    for name in df_wide_filtered.organism.unique():
        s, g, n, c = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c)
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sns.scatterplot(tmp, x='AF3 total residue plddt over 80 diff.', y='Foldseek bitscore diff.', alpha=0.3, ax=axes[0],
                    hue='Organism', legend=None, s=20)
    sns.scatterplot(tmp, x='AF3 total residue plddt over 80 diff.', y='IPR domain total length diff.', alpha=0.3, ax=axes[1],
                    hue='Organism', legend=None, s=20)
    sns.scatterplot(tmp, x='Foldseek bitscore diff.', y='IPR domain total length diff.', alpha=0.3, ax=axes[2], hue='Organism', s=20)
    for index in [0, 1, 2]:
        alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
        axes[index].annotate(
            alphabet[index],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3.0))
        axes[index].axvline(0, color='black', linestyle='-', lw=1)
        axes[index].axhline(0, color='black', linestyle='-', lw=1)

    sns.move_legend(axes[2], 'upper left', bbox_to_anchor=(1, 1), labels=labels, fontsize=10, title_fontsize=14)
    plt.suptitle('maxabs transformed')
    plt.tight_layout()
    plt.savefig('pairplot-best-scores-and-disorder-scatter-maxabs-diffs.png', dpi=300)
    plt.show()

def pairplot_diffs(df):
    df['diff_IPR_total'] = df["IPR_total_new"] - df["IPR_total_old"]
    sns.set_context('talk')
    tmp = df[['organism', 'diff_total_residue_plddt_over_80',
              'diff_bits', 'diff_IPR_total']].rename(columns={'organism': 'Organism',
                                                              'diff_total_residue_plddt_over_80': 'AF3 count of pLDDT above 80 diff.',
                                                              'diff_bits': 'Foldseek bitscore diff.',
                                                              'diff_IPR_total': 'IPR total diff.',
                                                              })
    labels = []
    for name in df_wide_filtered.organism.unique():
        s, g, n, c = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c)
    fig, axes = plt.subplots(1,3, figsize=(20,6))
    sns.scatterplot(tmp, x ='AF3 count of pLDDT above 80 diff.', y = 'Foldseek bitscore diff.' , alpha=0.3, ax=axes[0],
                    hue='Organism', legend=None, s=20)
    sns.scatterplot(tmp, x='AF3 count of pLDDT above 80 diff.', y='IPR total diff.', alpha=0.3, ax=axes[1],
                    hue='Organism', legend=None, s=20)
    sns.scatterplot(tmp, x='Foldseek bitscore diff.', y= 'IPR total diff.', alpha=0.3, ax=axes[2], hue='Organism', s=20)
    for index in [0,1,2]:
        alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
        axes[index].annotate(
            alphabet[index],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8,edgecolor='none', pad=3.0))
        axes[index].axvline(0, color='black', linestyle='-', lw=1)
        axes[index].axhline(0, color='black', linestyle='-', lw=1)

    sns.move_legend(axes[2], 'upper left', bbox_to_anchor=(1, 1), labels=labels,  fontsize=10, title_fontsize=14)
    plt.tight_layout()
    plt.savefig('pairplot-best-scores-and-disorder-scatter-diffs.png', dpi=300)
    plt.show()

pairplot_diffs(df_wide_filtered)
#pairplot_maxabs_diffs()
#pairplot_best(df_wide_filtered)
#af3_paired_disorder_plot(df_wide_filtered)