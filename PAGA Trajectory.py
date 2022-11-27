#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import fa2


# In[2]:


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/paul15.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')


# In[3]:


adata = sc.read_h5ad('tumor1_v2.h5ad')


# In[4]:


adata


# In[5]:


adata.X = adata.X.astype('float64')


# In[6]:


sc.pp.recipe_zheng17(adata)


# In[7]:


sc.tl.pca(adata, svd_solver='arpack')


# In[8]:


sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)


# In[9]:


sc.pl.draw_graph(adata, color='cell_type', legend_loc='on data')


# In[10]:


sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')


# In[11]:


sc.tl.draw_graph(adata)


# In[12]:


sc.pl.draw_graph(adata, color='cell_type', legend_loc='on data')


# In[13]:


sc.tl.leiden(adata, resolution=0.1)


# In[14]:


sc.tl.paga(adata, groups='leiden')


# In[16]:


sc.pl.paga(adata, color=['cell_type','leiden'])


# In[17]:


#skipped line 16 in tutorial:
#sc.pl.paga(adata, color=['louvain', 'Itga2b', 'Prss34', 'Cma1'])
adata.obs['leiden'].cat.categories


# In[21]:


#instructed to skip line 18
#skip line 19 too?:
adata.obs['leiden_anno'].cat.categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10/Ery', '11', '12',
       '13', '14', '15', '16/Stem', '17', '18', '19/Neu', '20/Mk', '21', '22/Baso', '23', '24/Mo']


# In[18]:


#skip 20 too?:
sc.tl.paga(adata, groups='leiden_anno')


# In[19]:


sc.pl.paga(adata, threshold=0.03, show=False)


# In[20]:


sc.tl.draw_graph(adata, init_pos='paga')


# In[21]:


sc.pl.draw_graph(adata, color=['cell_type','leiden'], legend_loc='on data')


# In[22]:


pl.figure(figsize=(8, 2))
for i in range(28):
    pl.scatter(i, 1, c=sc.pl.palettes.zeileis_28[i], s=200)
pl.show()


# In[33]:


zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['leiden_anno_colors'])
# error on second line of line 25
# original is 'louvain_anno_colors' but doesn't work either


# In[23]:


adata.uns['leiden_anno_colors'] = new_colors


# In[24]:


sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)


# In[ ]:




