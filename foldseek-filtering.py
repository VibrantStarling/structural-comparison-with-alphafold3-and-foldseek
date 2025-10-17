import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from sklearn.preprocessing import MaxAbsScaler


df_wide_filtered = pd.read_csv("E:/Dropbox/1-Veupath-files/PROJECTS/alphafold-other-tests/toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3.csv", sep = '\t',  low_memory=False)

thresholds = [1,0.05,0.01,0.001,0.0001,0.00001]
scores = ['lddt',  'bits', 'alntmscore']

masked_df = df_wide_filtered[['organism','lddt_new', 'bits_new', 'alntmscore_new', 'evalue_new',
                              'lddt_old', 'bits_old', 'alntmscore_old', 'evalue_old']]

masked_dict = {}
for threshold in thresholds:
    for name in scores:
        masked_dict[name + '_new_' + str(threshold)] = masked_df[name + '_new'].mask(masked_df['evalue_new'] > threshold, 0)
        masked_dict[name + '_old_' + str(threshold)] = masked_df[name + '_old'].mask(masked_df['evalue_old'] > threshold, 0)

tmp = pd.DataFrame(masked_dict)
tmp['Organism'] = masked_df['organism']

for threshold in thresholds:
    for v in scores:
        tmp.loc[:, 'diff_' + v+'_'+str(threshold)] = tmp.loc[:, v + '_new_'+str(threshold)] - tmp.loc[:, v + '_old_'+str(threshold)]
    for variable in scores:
        tmp['diff_'+variable+'_'+str(threshold)] = tmp['diff_'+variable+'_'+str(threshold)].mask(tmp['diff_'+variable+'_'+str(threshold)]>0, 1).mask(tmp['diff_'+variable+'_'+str(threshold)]<0, -1)

sns.set_style("white")
sns.set_context('talk')
fig, axs = plt.subplots(3,3, figsize=(15,20), sharey=False)
sns.set_palette(['#b5d1ae', '#568b87', '#122740'])

for row, score in enumerate(scores):
    scores_2 = ['diff_'+ score + '_' + str(threshold) for threshold in thresholds]
    scores_2.append('Organism')
    tmp2 = tmp[scores_2]
    # sort and plot data
    for index, organism in enumerate(list(tmp['Organism'].unique())):
        tmp3 = tmp2.loc[tmp['Organism'] == organism].melt(id_vars='Organism').groupby(
            ['variable', 'value']).size().reset_index().pivot(columns='variable', index='value', values=0)
        tmp3[['diff_'+score+'_1', 'diff_'+score+'_0.05', 'diff_'+score+'_0.01', 'diff_'+score+'_0.001',
              'diff_'+score+'_0.0001', 'diff_'+score+'_1e-05']].T.plot(kind='bar', stacked=True, ax=axs[row, index], width=1, legend=False)
        x_labels = list(tmp3.columns)
        axs[row, index].set_xlabel(None)
for index, organism in enumerate(list(tmp['Organism'].unique())):
        s, g, n, c = organism.split(' ')
        print(index)
        organism = "$\ " + s + " $" + " " + "$\ " + g + " $" + " " + n + " " + c
        axs[0, index].set_title(organism, pad=20)
for i, ax in enumerate(axs.flat):
    alphabet = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)', 'p)']
    ax.annotate(
        alphabet[i],
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', fontstyle='italic', verticalalignment='top',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3.0))
axs[0,2].legend(loc='upper left', bbox_to_anchor = (1,1))
plt.tight_layout()
plt.savefig('stacked-bar-thresholds-all.png', dpi=300)
plt.show()

