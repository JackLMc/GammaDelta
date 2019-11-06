
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import os
import glob
import pdb
import re
from collections import defaultdict, Counter, OrderedDict
import json
from itertools import imap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import subprocess
import pickle
from HTSeq import SAM_Reader
import Levenshtein
import pickle
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Levenshtein import distance
import pandas as pd
import seaborn as sns
import numpy as np
import scipy as sp
from scipy.stats import mannwhitneyu as mwu
from matplotlib import pyplot as plt
from matplotlib_venn import venn2_circles, venn2
import random
import networkx as nx
from pylab import rcParams
from IPython.display import Image
import copy
from operator import itemgetter
import json
sns.set_style("ticks")


# I want to do a first-pass look at clonotypes and expression from the data generated in Seurat on my local mac. This probably won't be the final way that I look at the data but it's quick and I want some things to show in the group meetings tomorrow and on Monday.

# In[2]:


tsne = pd.read_csv("../FirstData_Aug2016/tSNE_coords.txt",sep=" ")


# In[3]:


tsne.head()


# In[4]:


tsne.index = [x[1:] for x in tsne.index]


# In[5]:


tsne.head()


# In[6]:


tsne.index = ["_".join(x.split("_")[0:2]) + "#" + x.split("_")[2] for x in tsne.index]


# In[7]:


tsne.head()


# In[8]:


tpms = pd.read_csv("../FirstData_Aug2016/log_tpms_with_gene_names.txt", sep=" ")


# In[9]:


tpms.columns = [x[1:] for x in tpms.columns.values]


# In[10]:


tpms.columns = ["_".join(x.split("_")[0:2]) + "#" + x.split("_")[2] for x in tpms.columns.values]


# In[11]:


tpms.head()


# In[12]:


name_map = pd.read_csv("../FirstData_Aug2016/name_map.txt", sep="\t", index_col=0)


# In[13]:


name_map


# In[14]:


tsne = pd.merge(tsne, name_map.set_index('fastq_name'), left_index=True, right_index=True)


# In[15]:


name_map.query('fastq_name=="20427_6#98"')


# In[16]:


tsne['plate'] = ["_".join(x.split("_")[0:2]) for x in tsne['well']]


# In[17]:


tsne.head()


# In[18]:


def get_sort_class(well):
    well = well.split("_")[-1]
    number = int(well[1:])
    if number in [1,3,5,7,9,11]:
        return('naive')
    else:
        return('EMRA')


# In[19]:


tsne['sort_class'] = [get_sort_class(x) for x in tsne['well']]


# In[20]:


plt.rcParams['figure.figsize']=(10,10)


# In[21]:


sns.lmplot('tSNE_1', 'tSNE_2', data=tsne, hue='plate', fit_reg=False, size=10, scatter_kws={'s':50})


# In[22]:


sns.lmplot('tSNE_1', 'tSNE_2', data=tsne, hue='sort_class', fit_reg=False, size=10, scatter_kws={'s':50})


# In[23]:


clonotypes = []
#lanes_to_donors = {'8':'28', '7':'31', '6':'29'}
with open("../FirstData_Aug2016/clonotypes.txt") as f:
    for line in f:
        line = line.rstrip()
        line = line.split(",")
        clonotypes.append(line)
        


# In[24]:


for clonotype in clonotypes:
    if len(clonotype)>1:
        under_plot = tsne.drop(clonotype)
        clonotype_plot = tsne.ix[clonotype]
        sns.lmplot('tSNE_1', 'tSNE_2', data=under_plot, hue='sort_class', fit_reg=False, size=10, scatter_kws={'s':25})
        sns.regplot('tSNE_1', 'tSNE_2', data=clonotype_plot, fit_reg=False, scatter_kws={'s':50, 'c':'red'})
        print(Counter(clonotype_plot['sort_class']))
        plt.show()


# In[25]:


genes_of_interest = ['CD27', 'CD28', 'IL7R', 'CCR7', 'SELL', 'FCGR3A','GZMB', 'GZMA', 'PRF1', 'TCF7', 'CX3CR1', 'LEF1']


# In[26]:


tpms.head()


# In[27]:


tpms.ix['IL1B']


# In[28]:


cytokine_list = [['IL2'], ['IL4'], ['IL5'], ['IL6'], ['IL9'], ['CXCL8'], ['IL10'], ['IL12A', 'IL12B'], ['IL12B'], ['IL17A'], ['IL17F'], ['IL18'], 
                 ['IL23A', 'IL12B'], ['IL27'], ['IL33'], ['CFS2'], ['TNF'], ['IFNG'], 
                 ['IFNA1', 'IFNA2', 'IFNA4', 'IFNA5','IFNA6', 'IFNA7', 'IFNA8','IFNA10', 'IFNA13', 'IFNA14', 'IFNA16',
                  'IFNA17', 'IFNA21'], ['FAS'], ['FASLG'], ['CCL2'], ['IL1B']]


# In[29]:


tsne.head()


# In[32]:


sort_details = tsne[:,2:]


# In[ ]:


for c in cytokine_list:
    df = tsne.copy()
    t = tpms.ix[c]
    
    t = pd.Series(t.sum(0), name="_".join(c))
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    plt.scatter(df['tSNE_1'], df['tSNE_2'], c=df["_".join(c)], cmap='YlOrRd')
    print ("_".join(c))
    plt.colorbar()
    plt.show()


# In[ ]:


plt.scatter(df['tSNE_1'], df['tSNE_2'], c=df["_".join(c)], cmap='YlOrRd')
plt.title="Test"


# In[ ]:


interesting_cytokine_list = [['IL2'], ['IL4'], ['IL5'], ['IL6'], ['CXCL8'], ['IL10'], ['IL18'], 
                 ['IL23A', 'IL12B'], ['IL33'], ['TNF'], ['IFNG'], 
                 ['IFNA1', 'IFNA2', 'IFNA4', 'IFNA5','IFNA6', 'IFNA7', 'IFNA8','IFNA10', 'IFNA13', 'IFNA14', 'IFNA16',
                  'IFNA17', 'IFNA21'], ['FAS'], ['FASLG'], ['IL1B']]


# In[ ]:


len(interesting_cytokine_list)


# In[ ]:


plt.rcParams['figure.figsize']=(25,15)


# In[ ]:


name


# In[ ]:


df.shape


# In[ ]:


c


# In[ ]:


f, axs = plt.subplots(ncols=5, nrows=3)
for i in range(len(interesting_cytokine_list)):
    row = i / 5
    col = i % 5
    ax = axs[row, col]
    df = tsne.copy()
    c = interesting_cytokine_list[i]
    t = tpms.ix[c]
    name = "_".join(c)
    t = pd.Series(t.sum(0), name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    vmin = min(df[name])
    vmax = max(df[name])
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("cytokines.pdf")


# In[ ]:


nc = 5
nr = 3
f, axs = plt.subplots(ncols=nc, nrows=nr)
for i in range(len(interesting_cytokine_list)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = interesting_cytokine_list[i]
    t = tpms.ix[c]
    name = "_".join(c)
    t = pd.Series(t.sum(0), name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    vmin = 0
    vmax = 7.5
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("cytokines_same_scale.pdf")


# In[ ]:


exhaustion_markers = ['PDCD1', 'LAG3', 'HAVCR2', 'IL2RB', 'KLRG1']
tpms.ix[exhaustion_markers]


# In[ ]:


genes = exhaustion_markers
nc = 3
nr = 2
f, axs = plt.subplots(ncols=nc, nrows=nr)
for i in range(len(genes)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = genes[i]
    t = tpms.ix[c]
    name = c
    t = pd.Series(t, name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    vmin = min(df[name])
    vmax = max(df[name])
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("exhaustion_markers.pdf")


# In[ ]:


genes = exhaustion_markers
nc = 3
nr = 2
f, axs = plt.subplots(ncols=nc, nrows=nr)
vmin = min(tpms.ix[genes].min())
vmax = max(tpms.ix[genes].max())
for i in range(len(genes)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = genes[i]
    t = tpms.ix[c]
    name = c
    t = pd.Series(t, name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("exhaustion_markers_same_scale.pdf")


# In[ ]:


TFs = ['EOMES', 'TBX21', 'PRDM1', 'RORC', 'BCL6', 'GATA3']
tpms.ix[TFs]


# In[ ]:


genes = TFs
nc = 3
nr = 2
f, axs = plt.subplots(ncols=nc, nrows=nr)
for i in range(len(genes)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = genes[i]
    t = tpms.ix[c]
    name = c
    t = pd.Series(t, name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    vmin = min(df[name])
    vmax = max(df[name])
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("TFs.pdf")


# In[ ]:


genes = TFs
nc = 3
nr = 2
f, axs = plt.subplots(ncols=nc, nrows=nr)
vmin = min(tpms.ix[genes].min())
vmax = max(tpms.ix[genes].max())
for i in range(len(genes)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = genes[i]
    t = tpms.ix[c]
    name = c
    t = pd.Series(t, name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("TFs_same_scale.pdf")


# In[ ]:


activation_markers = ['IL18R1', 'SELL', 'IL2RB', 'CXCR3']
tpms.ix[activation_markers]


# In[ ]:


genes = activation_markers
nc = 2
nr = 2
f, axs = plt.subplots(ncols=nc, nrows=nr)
for i in range(len(genes)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = genes[i]
    t = tpms.ix[c]
    name = c
    t = pd.Series(t, name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    vmin = min(df[name])
    vmax = max(df[name])
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("activation_markers.pdf")


# In[ ]:


genes = activation_markers
nc = 2
nr = 2
f, axs = plt.subplots(ncols=nc, nrows=nr)
vmin = min(tpms.ix[genes].min())
vmax = max(tpms.ix[genes].max())
for i in range(len(genes)):
    row = i / nc
    col = i % nc
    ax = axs[row, col]
    df = tsne.copy()
    c = genes[i]
    t = tpms.ix[c]
    name = c
    t = pd.Series(t, name=name)
    #t = np.exp(t) - 1
    df = pd.concat([df, t],1)
    positive = df.query("{}>0".format(name))
    negative = df.drop(positive.index.values)
    ax.scatter(negative['tSNE_1'], negative['tSNE_2'], c=negative[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.scatter(positive['tSNE_1'], positive['tSNE_2'], c=positive[name], cmap='YlOrRd', s=125, vmin=vmin, vmax=vmax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if 'IFNA1' in name:
        name = 'IFNA'
    ax.set_title(name)
plt.savefig("activation_markers_same_scale.pdf")

