import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import scipy
import math
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib import *
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42
ad = sc.read('./whole_cell_pfc.h5ad')
ad.X = ad.layers['counts']
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad,base=2)
exp_data = ad.to_df()
exp_data['age_binary'] = ad.obs['age_binary']#age_binary contains two values 'young' or 'old'
plot_data = exp_data.groupby('age_binary').apply(lambda x:x.mean(0))
neuron_dev_genes = pd.read_table('../neuron_develop_genes.txt',header=None)
neuron_dev_genes = neuron_dev_genes[0].tolist()
plot_data = plot_data.loc[:,plot_data.columns.isin(neuron_dev_genes)]
# Function to find distance
def shortest_distance(x1, y1, a=-1, b=1, c=0):    
    d = abs((a * x1 + b * y1 + c)) / (math.sqrt(a * a + b * b))
    return d
for i in np.unique(ad.obs['label']):
    ad_temp = ad[ad.obs['label']==i]
    exp_data = ad_temp.to_df()
    exp_data['age_binary'] = ad_temp.obs['age_binary']
    plot_data = exp_data.groupby('age_binary').apply(lambda x:x.mean(0))
    plot_data = plot_data.loc[:,plot_data.columns.isin(neuron_dev_genes)]
    r_value,p_value = scipy.stats.pearsonr(plot_data.loc['Young',],plot_data.loc['Old',])
    distance_series = pd.Series()
    for j in plot_data.columns:
        d= shortest_distance(x1=plot_data.loc['Young',j],y1=plot_data.loc['Old',j])
        distance_series[j]=d
    distance_series = distance_series.sort_values(ascending=False)
    sele_gene = distance_series.index[0:20].tolist()
    fig = plt.figure(figsize=(4,4))
    ax = plt.axes()
    red_gene=[]
    for index in plot_data.columns:
        if plot_data.loc['Old',index] - plot_data.loc['Young',index] > 1:
            red_gene.append(index)
        else:
            if plot_data.loc['Old',index] - plot_data.loc['Young',index] < -1:
                red_gene.append(index)
    for gene_index in red_gene:
        if gene_index in sele_gene:
            sele_gene.remove(gene_index)
    blue_gene = sele_gene[0:5]
    color_list=[]
    for gene_index in plot_data.columns:
        if gene_index in red_gene:
            color_list.append('red')
        else:
            if gene_index in blue_gene:
                color_list.append('blue')
            else:
                color_list.append('grey')
    plt.scatter(plot_data.loc['Young',],plot_data.loc['Old',],s=1,color=color_list)
    # add label
    plot_text_gene = red_gene + blue_gene
    texts = [plt.text(plot_data.loc['Young',gene_index],plot_data.loc['Old',gene_index],gene_index) for gene_index in plot_text_gene]
    plt.text(0.3,3.5,f'r={round(r_value,3)}')
    plt.text(0.3,3,f'p={round(p_value,3)}')
    adjust_text(texts)
    plt.xlim([0,4.5])
    plt.ylim([0,4.5])
    plt.axline((0, 0), (4, 4),color='black')
    plt.axline((0,1),(3,4),linestyle='--',color='black',linewidth=0.5)
    plt.axline((1,0),(4,3),linestyle ='--',color='black',linewidth=0.5)
    ax.set_xticks([0,1,2,3,4])
    ax.set_yticks([0,1,2,3,4])
    plt.xlabel('Young',size=12)
    plt.ylabel('Old',size=12)
    plt.title(i)
    plt.show()
    fig.savefig(f'./Figures/PFC_{i}_compare_young_old_gene.pdf')
