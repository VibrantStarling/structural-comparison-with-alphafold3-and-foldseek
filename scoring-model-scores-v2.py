import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from sklearn.preprocessing import MaxAbsScaler

#df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED.tsv", sep = '\t', low_memory=False)
df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3.csv", sep = '\t',  low_memory=False)

print(df_wide_filtered)

evalue_is_sig = 0.05

increase_is_good = [
                    'diff_ptm',
                    'diff_residue_plddt_mean',
                    'diff_total_residue_plddt_over_80',
                    'diff_bits',
                    # 'diff_fident',
                    'diff_lddt',
                    'diff_alntmscore',
                    'diff_Foldseek_transformed_evalue',
                    #'diff_PSAURON_in_frame_score',
                    # 'diff_IPR_transformed_min_evalue',
                    'diff_IPR_total',
                    'diff_IPR_max'
                    ]
decrease_is_good = ['diff_evalue',
                    'diff_pae_mean',
                    'diff_mean_disorder_metapredict3',
                    'diff_IPR_min_evalue'
                    ]
other = ['organism']


PTM = {-1 : (0,0.5),
       1 : (0.5,1)}


phLDDT = {-1:[0,0.5],
      1:[0.5,1]}


psauron = {-1:[0,0.7],0:[0.7,0.95],
      1:[0.95,1]}


print(df_wide_filtered['evalue_old'].max())

# set Foldseek scores to 0 where evalue is >0.001, and the evalue to -1
values = ['lddt', 'fident', 'bits', 'alntmscore', 'Foldseek_transformed_evalue']

for name in values:
    df_wide_filtered.loc[df_wide_filtered['evalue_new'] > evalue_is_sig, name+'_new'] = 0
    df_wide_filtered.loc[df_wide_filtered['evalue_old'] > evalue_is_sig, name+'_old'] = 0
    df_wide_filtered[name+'_new'].fillna(0, inplace=True)
    df_wide_filtered[name+'_old'].fillna(0, inplace=True)

values = ['IPR_total', 'IPR_max']
for name in values:
    df_wide_filtered[name+'_new'].fillna(0, inplace=True)
    df_wide_filtered[name+'_old'].fillna(0, inplace=True)


# sort out null values in evalue
df_wide_filtered.loc[df_wide_filtered['evalue_new'] > evalue_is_sig, 'evalue_new'] = 10
df_wide_filtered.loc[df_wide_filtered['evalue_old'] > evalue_is_sig, 'evalue_old'] = 10
df_wide_filtered['evalue_new'].fillna(10, inplace=True)
df_wide_filtered['evalue_old'].fillna(10, inplace=True)
df_wide_filtered.loc[df_wide_filtered['evalue_new'] == 0, 'evalue_new'] = df_wide_filtered['evalue_new'].min()
df_wide_filtered.loc[df_wide_filtered['evalue_old'] == 0, 'evalue_old'] = df_wide_filtered['evalue_new'].min()



# recalculate negative logs
df_wide_filtered['Foldseek_transformed_evalue_new'] = np.log10(df_wide_filtered['evalue_new'])*-10
df_wide_filtered['Foldseek_transformed_evalue_old'] = np.log10(df_wide_filtered['evalue_old'])*-10
#df_wide_filtered['diff_Foldseek_transformed_evalue'] = df_wide_filtered['Foldseek_transformed_evalue_new'] - df_wide_filtered['Foldseek_transformed_evalue_old']



values = ['lddt', 'fident', 'bits', 'alntmscore', 'Foldseek_transformed_evalue', 'evalue',
          'IPR_evalue_min', 'IPR_transformed_evalue_min','IPR_total','IPR_max','PSAURON_in_frame_score','pae_mean','ptm',
          'residue_plddt_median','total_residue_plddt_over_80']

print(df_wide_filtered.columns)
# recalculate diffs
for v in values:
    df_wide_filtered.loc[:,'diff_'+v] = df_wide_filtered.loc[:, v+'_new'] - df_wide_filtered.loc[:, v+'_old']

#invert decrease_is_good = ['diff_evalue', 'diff_pae_mean']
df_wide_filtered[decrease_is_good] = df_wide_filtered[decrease_is_good] * -1


# redefine decrease is good to remove raw evalues
decrease_is_good = ['diff_pae_mean',
                    'diff_mean_disorder_metapredict3',
                    ]
df_wide_filtered.drop(columns=['diff_evalue',
                    'diff_IPR_min_evalue'], inplace=True)

# renaming everything for publication
df_wide_filtered = df_wide_filtered.rename(columns={
                    'organism':'Organism',
                    'diff_ptm':'AF3 pTM diff.',
                    'diff_residue_plddt_median':'AF3 plDDT median diff.',
                    'diff_total_residue_plddt_over_80':'AF3 total plDDT over 80 diff.',
                    'diff_bits':'Foldseek bitscore diff.',
                    'diff_lddt':'Foldseek LDDT diff.',
                    'diff_alntmscore':'Foldseek alntmscore diff.',
                    'diff_Foldseek_transformed_evalue':'Foldseek log10(evalue)*-10 diff.',
                    'diff_IPR_total':'IPR domain total length diff.',
                    'diff_IPR_max':'IPR domain max. length diff.',
                    'diff_pae_mean':'AF3 mean PAE diff.',
                    'diff_mean_disorder_metapredict3':'Metapredict3 mean diff.'
})

increase_is_good = [
                    'AF3 pTM diff.',
                    'AF3 plDDT median diff.',
                    'AF3 total plDDT over 80 diff.',
                    'Foldseek bitscore diff.',
                    'Foldseek LDDT diff.',
                    'Foldseek alntmscore diff.',
                    'Foldseek log10(evalue)*-10 diff.',
                    'IPR domain total length diff.',
                    'IPR domain max. length diff.'
                    ]
decrease_is_good = [#'Metapredict3 mean disorder diff.',
                    'AF3 mean PAE diff.'
                    ]
other = ['Organism']


def get_all():
    # apply MaxAbsScaler
    test = df_wide_filtered[other + decrease_is_good +increase_is_good]
    X = df_wide_filtered[decrease_is_good + increase_is_good ]
    y = df_wide_filtered[other]
    print(X)
    transformer = MaxAbsScaler().fit(X)
    X_t = transformer.transform(X)
    df_transformed = pd.DataFrame(X_t, columns=decrease_is_good + increase_is_good )
    df_transformed['Organism'] = test['Organism']

    # apply mask to differences
    for name in decrease_is_good+ increase_is_good:
        test[name] = test[name].mask(test[name]>0, 1).mask(test[name]<0, -1)
    print(test)
    return test, df_transformed

def get_some(list_of_variables):
    # apply MaxAbsScaler
    test = df_wide_filtered[['Organism'] + list_of_variables]
    x = df_wide_filtered[list_of_variables]
    print(x)
    transformer = MaxAbsScaler().fit(x)
    x_t = transformer.transform(x)
    df_transformed = pd.DataFrame(x_t, columns=list_of_variables)
    df_transformed['Organism'] = test['Organism']
    # apply mask to differences
    for variable in list_of_variables:
        test[variable] = test[variable].mask(test[variable]>0, 1).mask(test[variable]<0, -1)
    print(test)
    return test, df_transformed

test, df_transformed = get_all()


def generate_win_lose_sum_cols():
    test, df_transformed = get_all()
    test_3, df_transformed_3 = get_some(['AF3 total plDDT over 80 diff.','Foldseek bitscore diff.', 'IPR domain total length diff.'])
    test_AF3Fold, df_transformed_AF3Fold = get_some(['AF3 total plDDT over 80 diff.','Foldseek bitscore diff.'])
    test_AF3IPR, df_transformed_AF3IPR = get_some(['AF3 total plDDT over 80 diff.', 'IPR domain total length diff.'])

    df_transformed['Sum of all'] = df_transformed[increase_is_good + decrease_is_good].sum(axis=1)
    df_transformed['Sum of AF3+Foldseek+IPR'] = df_transformed_3.drop(columns='Organism').sum(axis=1)
    df_transformed['Sum of AF3+Foldseek'] = df_transformed_AF3Fold.drop(columns='Organism').sum(axis=1)
    df_transformed['Sum of AF3+IPR'] = df_transformed_AF3IPR.drop(columns='Organism').sum(axis=1)

    for variable in ['Sum of all', 'Sum of AF3+Foldseek+IPR',
                     'Sum of AF3+Foldseek', 'Sum of AF3+IPR']:
        test[variable] = df_transformed[variable].mask(df_transformed[variable]>0, 1).mask(df_transformed[variable]<0, -1)
    return test, df_transformed

test, df_transformed = generate_win_lose_sum_cols()

def plot_ecdf():
    print(df_transformed.melt(id_vars='Organism'))
    labels = []
    for name in df_transformed.Organism.unique():
        s, g , n, c  = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c)
    print(labels)
    g = sns.displot(df_transformed.melt(id_vars='Organism'), x='value', col='variable', col_wrap=3, hue='Organism',
                hue_order=df_transformed.Organism.unique(), kind="ecdf", facet_kws=dict(legend_out=False))
    sns.move_legend(g,  "upper left", bbox_to_anchor=(0.8, 0.95), labels=labels)
    plt.show()
    return

#plot_ecdf()

def plot_heatmap(df):
    cmap = sns.diverging_palette(100, 5, s=100, l=40, center="dark", n=9)
    fig, ax = plt.subplots(1,3)
    for index, organism in enumerate(list(df['Organism'].unique())):
        tmp = df.loc[df['Organism'] == organism, increase_is_good + decrease_is_good].fillna(0)
        tmp.columns = [i.capitalize().replace('_',' ') for i in list(tmp.columns)]
        s, g, n, c = organism.split(' ')
        label = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        sns.heatmap(tmp,
                    cmap=cmap, ax=ax[index], yticklabels=False, xticklabels=True, center=0,
                    vmin=-1, vmax=1)
        ax[index].set_title(label)
    plt.show()
    return

#plot_heatmap(df)


def plot_cluster(df, name, organism):
    cmap = sns.diverging_palette(240, 10, s=70, l=60, center="dark", n=100)
    s, g , n, c  = organism.split(' ')
    o = "$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c
    print(o)
    tmp = df[decrease_is_good + increase_is_good].fillna(0)
    tmp = tmp.rename(columns={'diff_ptm':'diff_pTM'})
    cg = sns.clustermap(tmp, row_cluster=True, col_cluster=False,  yticklabels=False, xticklabels=True, center=0,
                        cmap=cmap, vmin=-1, vmax=1, dendrogram_ratio=(.1, .2), tree_kws={"linewidths": 0.},
                        cbar_pos=(0.02, .3, .03, .5))
    plt.suptitle(o + ', ' + name.lower(), x=.55, y = 0.9)
    plt.savefig(name+'_'+s+g+'.png', dpi=300)
    plt.savefig(name + '_' + s + g + '.svg', format='svg',dpi=300)
    plt.show()
    return


def plot_those_clusters():
    #for organism in list(test['Organism'].unique()):
     #   plot_cluster(test.loc[test['Organism']==organism], 'Masked', organism)

    for organism in list(df_transformed['Organism'].unique()):
        plot_cluster(df_transformed.loc[df_transformed['Organism']==organism], 'MaxAbs Transformed', organism)

plot_those_clusters()

def plot_hist():
    labels = []
    tmp = df_transformed[['Organism','Sum of AF3+Foldseek+IPR',
                     'Sum of AF3+Foldseek', 'Sum of AF3+IPR']].melt(id_vars='Organism')
    #tmp.rename(columns={'value':'Sum of differences'})
    print(tmp)
    sns.set_style('white')
    for name in tmp.Organism.unique():
        s, g , n, c  = name.split(' ')
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c)
    g = sns.displot(tmp, row='variable', x='value',
                    col='Organism', hue='Organism',facet_kws=dict(legend_out=False, sharex=True, sharey=False),
                    legend=None, fill=True, palette='colorblind', height=4, aspect=1.2)
    g.set_titles(template="")
    alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
    y_labels = ['AF3+Foldseek+IPR, count', 'AF3+Foldseek, count', 'AF3+IPR, count']
    index = 0
    for ax in g.axes.flat:
        ax.annotate(
            alphabet[index],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize=16, fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
        index += 1
    for i, ax in enumerate(g.axes[:, 0]):
        ax.set_ylabel(y_labels[i], fontsize=16)
    for i, ax in enumerate(g.axes[2, :]):
        ax.set_xlabel('Sum of differences', fontsize=16)
    def specs(x, **kwargs):
        #plt.axvline(x.mean(), color='black', ls='-', lw=1)
        plt.axvline(0, color='black', linestyle='--', lw=1)
    g.map(specs, 'value')
    for index, organism in enumerate(labels):
        g.axes[0, index].set_title(organism, fontdict = {'fontsize': 16})
    plt.tight_layout()
    plt.savefig("sum-of-diffs-displot.png", dpi=300)
    plt.show()
    return


def plot_hist_2():
    print(test.melt(id_vars='Organism'))
    test['Sum of masked differences'] = test[decrease_is_good + increase_is_good].sum(axis=1)
    labels = []
    for name in df_transformed.Organism.unique():
        s, g , n, c  = name.split(' ')
        print(n)
        labels.append("$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c)
    g = sns.displot(test[['Organism','Sum of masked differences']], x='Sum of masked differences', col='Organism',
                col_wrap=3, hue='Organism', kind="hist")
    sns.move_legend(g, "upper left", bbox_to_anchor=(0.8, 0.95), labels=labels)
    plt.show()
    return

plot_hist()
#plot_hist_2()


df_transformed.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-2-CHANGED-FILTERED-3-maxabs.tsv", sep = '\t', index=None)
test.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-FILTERED-3-win-lose-masked.tsv", sep = '\t', index=None)
df_wide_filtered.to_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-FILTERED-3-evalue-masked.tsv", sep = '\t', index=None)


#sns.set_context('paper')

def plot_pairs():
    p = sns.pairplot(df_transformed, hue='Organism')
    labels = []
    for name in df_transformed.Organism.unique():
            s, g , n, c  = name.split(' ')
            labels.append("$\ " + s + " $" + " " + "$\ " + g + " $"+ " " + n + " " + c)
    sns.move_legend(p, "upper left", bbox_to_anchor=(0.9, 0.9), labels=labels)
    plt.savefig('differences_maxabs_pairplot.png', dpi=300)
    plt.show()
    plt.clf()



## plot stacked bar
def plot_stacked():
    sns.set_style("white")
    fig, ax = plt.subplots(1,3, figsize=(15,8), sharey=False)
    sns.set_palette(['#b5d1ae', '#568b87', '#122740'])
    sns.set_style('dark')
    for index, organism in enumerate(test.Organism.unique()):
        tmp = test.loc[test['Organism']==organism].melt(id_vars = 'Organism').groupby(['variable','value']).size().reset_index().pivot(columns='variable', index='value', values=0)
        tmp.insert(4, ' ', value=np.nan)
        tmp.insert(9, '   ', value=np.nan)
        tmp.insert(12, '  ', value = np.nan)
        print(tmp)
        tmp.T.plot(kind='bar', stacked=True, ax = ax[index],  width=1, legend=False)
        x_labels = list(tmp.columns)
        print(x_labels)
        s, g, n, c = organism.split(' ')
        organism = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        ax[index].set_title(organism)
        alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
        ax[index].annotate(
            alphabet[index],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='large', fontstyle='italic', verticalalignment='top',
            bbox=dict(facecolor='none', edgecolor='none', pad=3.0))
        ax[index].set_xlabel(None)
    ax[2].legend(loc='upper left', bbox_to_anchor = (1,1))
    plt.tight_layout()
    plt.savefig('scored_stacked_bar_blue'+str(evalue_is_sig)+'.png', format='png', dpi=300)
    plt.show()

plot_stacked()

