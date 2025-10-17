import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3.csv", sep = '\t', low_memory=False)

sns.set(style='white', palette='magma')

organisms = list(df_wide_filtered.organism.unique())

print(df_wide_filtered.groupby('organism').count())

def plot_diffs():
    sns.set_context('talk')
    fig, ax = plt.subplots(1, 3, figsize=(15, 10), sharey=True)
    for index, organism in enumerate(organisms):
        tmp = df_wide_filtered.loc[df_wide_filtered['organism']==organism, ['diff_residue_plddt_50',
        'diff_residue_plddt_60',
        'diff_residue_plddt_70',
        'diff_residue_plddt_80',
        'diff_residue_plddt_90',]]
        tmp.columns = ['diff. ≥50 to <60',
        'diff. ≥60 to <70',
        'diff. ≥70 to <80',
        'diff. ≥80 to <90',
        'diff. ≥90']
        sns.boxplot(tmp, ax=ax[index])
        ax[index].set_xticks(ax[index].get_xticks(), ax[index].get_xticklabels(), rotation=90)
        s, g, n, c = organism.split(' ')
        o = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        ax[index].set_title(o)
        alphabet = ['a)', 'b)', 'c)']
        ax[index].annotate(
            alphabet[index],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top', fontfamily='sanserif',
            bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
    ax[0].set_ylabel('Difference in pLDDT count per protein pair')
    plt.tight_layout()
    plt.savefig('plddt-boxplots.png')
    plt.show()
    return


#plot_diffs()

def plot_diffs_stacked():
    sns.set_context('talk')
    #sns.set_style('dark')
    sns.set_palette(['#b5d1ae', '#568b87', '#122740'])
    fig, ax = plt.subplots(1, 3, figsize=(15, 8), sharey=False)
    for index, organism in enumerate(organisms):
        tmp = df_wide_filtered.loc[df_wide_filtered['organism']==organism, ['organism', 'diff_above_90','diff_above_80',
                                                                                       'diff_above_70','diff_above_60',
                                                                                       'diff_above_50']]

        for variable in ['diff_above_50','diff_above_60','diff_above_70','diff_above_80','diff_above_90']:
            tmp[variable] = tmp[variable].mask(tmp[variable] > 0, 1).mask(
                tmp[variable] < 0, -1)
        counts_df = tmp[['diff_above_50','diff_above_60','diff_above_70',
                         'diff_above_80','diff_above_90']].melt().groupby(
            by=['variable'], as_index=False).value_counts().pivot(
            index='value', columns='variable', values='count')
        print(counts_df)
        counts_df.T.plot(kind='bar', stacked=True, ax = ax[index],  width=1, legend=False)
        s, g, n, c = organism.split(' ')
        organism = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        ax[index].set_title(organism)
        alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
        ax[index].annotate(
            alphabet[index],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='1', edgecolor='none', pad=5.0, alpha=0.7))
        ax[index].set_xlabel(None)
    ax[2].legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig('plddt-stacked-barplots.png')
    plt.show()
    return


plot_diffs_stacked()

def plot_pointplot(df):
    df_melt = df_wide_filtered[['organism','residue_plddt_50_new', 'residue_plddt_50_old',
    'residue_plddt_60_new',	'residue_plddt_60_old',
    'residue_plddt_70_new',	'residue_plddt_70_old',
    'residue_plddt_80_new',	'residue_plddt_80_old',
    'residue_plddt_90_new',	'residue_plddt_90_old']].melt(id_vars='organism')

    df_melt.loc[df_melt['variable'].str.contains('new'), 'source'] = 'new'
    df_melt.loc[df_melt['variable'].str.contains('old'), 'source'] = 'old'
    df_melt['variable'] = df_melt['variable'].str.replace('_new','').str.replace('_old','')

    sns.catplot(df_melt, x='source', y='value', hue='variable', order=['old','new'], col='organism', kind='point')
    plt.show()
    return


#plot_pointplot(df_wide_filtered)


def plot_cumsum_1():
    fig, ax = plt.subplots(1,2,figsize=(12,8))
    tmp_new = df_wide_filtered[['residue_plddt_50_new',
        'residue_plddt_60_new',
        'residue_plddt_70_new',
        'residue_plddt_80_new',
        'residue_plddt_90_new']]
    tmp_new.columns = ['≥50 to <60',
        '≥60 to <70',
        '≥70 to <80',
        '≥80 to <90',
        '≥90']
    tmp_new.cumsum(axis=1).plot(kind='box', ax=ax[0])

    tmp_old = df_wide_filtered[['residue_plddt_50_old',
        'residue_plddt_60_old',
        'residue_plddt_70_old',
        'residue_plddt_80_old',
        'residue_plddt_90_old']]
    tmp_old.columns = ['≥50 to <60',
        '≥60 to <70',
        '≥70 to <80',
        '≥80 to <90',
        '≥90']
    tmp_old.cumsum(axis=1).plot(kind='box',ax=ax[1])

    ax[0].set_title('New models')
    ax[0].set_xticks(ax[0].get_xticks(), ax[0].get_xticklabels(), rotation=90)
    ax[1].set_title('Old models')
    ax[1].set_xticks(ax[1].get_xticks(), ax[1].get_xticklabels(), rotation=90)
    plt.show()
    return


#plot_cumsum_1()


def plot_cumsum_2():
    fig, ax = plt.subplots(1,2,figsize=(12,8))
    tmp_old = df_wide_filtered[['residue_plddt_50_old',
            'residue_plddt_60_old',
            'residue_plddt_70_old',
            'residue_plddt_80_old',
            'residue_plddt_90_old']]
    tmp_old.columns = ['≥50 to <60',
                       '≥60 to <70',
                       '≥70 to <80',
                       '≥80 to <90',
                       '≥90']
    tmp_new = df_wide_filtered[['residue_plddt_50_new',
            'residue_plddt_60_new',
            'residue_plddt_70_new',
            'residue_plddt_80_new',
            'residue_plddt_90_new']]
    tmp_new.columns = ['≥50 to <60',
                       '≥60 to <70',
                       '≥70 to <80',
                       '≥80 to <90',
                       '≥90']
    sns.ecdfplot(tmp_new.iloc[:, ::-1].cumsum(axis=1),
                  ax=ax[0], stat='count', legend = None)
    sns.ecdfplot(tmp_old.iloc[:, ::-1].cumsum(axis=1),
                  ax=ax[1], stat='count')
    sns.move_legend(ax[1], 'upper left', bbox_to_anchor = (1,1))
    ax[0].set_xlabel('Number of models')
    ax[0].set_xlim(0,2000)
    ax[0].set_title('New models')
    ax[1].set_xlabel('Number of models')
    ax[1].set_xlim(0, 2000)
    ax[1].set_title('Old models')
    plt.suptitle('Cumulative sums of pLDDT counts')
    plt.show()

#plot_cumsum_2()


def alt_thresholds():
    df_residues_old = df_wide_filtered[['residue_plddt_50_old',
            'residue_plddt_60_old',
            'residue_plddt_70_old',
            'residue_plddt_80_old',
            'residue_plddt_90_old']].iloc[:, ::-1].cumsum(axis=1)

    df_residues_old.columns = ['above_90_old','above_80_old','above_70_old','above_60_old','above_50_old']
    df_residues_new = df_wide_filtered[['residue_plddt_50_new',
            'residue_plddt_60_new',
            'residue_plddt_70_new',
            'residue_plddt_80_new',
            'residue_plddt_90_new']].iloc[:, ::-1].cumsum(axis=1)
    df_residues_new.columns = ['above_90_new','above_80_new','above_70_new','above_60_new','above_50_new']
    df_residues = df_residues_old.join(df_residues_new)

    fig, ax = plt.subplots(1,5,  sharey=True, sharex=True)
    for i, n in enumerate(['50','60','70','80','90']):
        sns.ecdfplot(df_residues[['above_'+n+'_new','above_'+n+'_old']], ax=ax[i], palette='tab10')
        sns.move_legend(ax[i], 'lower right')
        #ax[0].xlabel('')
    plt.show()

    df_diffs = pd.DataFrame(np.array(df_residues_new)-np.array(df_residues_old), columns = ['diff_above_90','diff_above_80',
                                                                                       'diff_above_70','diff_above_60',
                                                                                       'diff_above_50']).join(
        df_wide_filtered[['diff_bits', 'organism']].reset_index(drop=True))
    fig, ax = plt.subplots(1, 5, sharey=True, sharex=True)
    for i, n in enumerate(['50', '60', '70', '80', '90']):
        sns.scatterplot(df_diffs[['diff_above_' + n , 'diff_bits','organism']], x= 'diff_above_' + n , y =  'diff_bits',
                        hue='organism', palette='tab10', ax=ax[i])
    ax[0].get_legend().remove()
    ax[1].get_legend().remove()
    ax[2].get_legend().remove()
    ax[3].get_legend().remove()
    sns.move_legend(ax[4], 'upper left', bbox_to_anchor=(1,1))
    plt.show()
    return df_residues.join(df_diffs).drop(columns=['organism','diff_bits'])

#alt_thresholds()

#df_diffs = alt_thresholds()
#df_wide_filtered = df_wide_filtered.reset_index(drop=True).join(df_diffs.reset_index(drop=True))
#df_wide_filtered.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-2.csv",sep='\t', index=False)




