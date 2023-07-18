import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
#import scvelo as scv
import scanpy.external as sce
import os,re

import seaborn as sns
#from statannot import add_stat_annotation
import matplotlib.ticker as mtick
os.chdir('/share/pub/zhangyr/projects/cooperation/SJTU/ESCC')

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')


def preprocessing(adata,filter_cells, n_genes_by_counts_max = 60000, n_genes_by_counts_min = 200, pct_counts_mt_min = 10):
    sc.pp.filter_genes(adata, min_cells=0.001*len(adata.obs))
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['rb'] = adata.var_names.str.contains('^RP[SL]')
    adata.var['hb'] = adata.var_names.str.contains('^HB[APS]')
    adata.var['hsp'] = adata.var_names.str.contains('^HSP')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rb','hb','hsp'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    adata = adata[adata.obs.n_genes_by_counts > 200, :]
    #adata = adata[adata.obs.total_counts < 20000, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    pre_cell = len(adata.obs_names)
    adata = adata[~adata.obs.index.isin(filter_cells)]
    after_cell = len(adata.obs_names)
    print("{0} cells filtered by double and mix".format(pre_cell-after_cell) )
    return adata

def normalization(adata, n_top = 3000):
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=n_top)
    #adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'n_genes_by_counts', 
                      'pct_counts_mt'])
    sc.tl.pca(adata, svd_solver='arpack')
    return adata


def getTypeCells(adata, meta, cell_type):
    cells = list(set(meta[meta['leiden_name_old'] == cell_type].index) & set(adata.obs.index))
    res = adata[adata.obs.index.isin(cells)]
    return res

## 想强制删掉的细胞
filter_cells = []

for data in datalist:
    ## 单个样本预处理
    data = preprocessing(data,[])

adata_con = datalist[0].concatenate([1:length(datalist)])

## 多个样本合起来标准化
adata_con = normalization(adata_con,min_cells=3)
## 多个样本去批次效应
sc.external.pp.bbknn(adata_con, batch_key='batch_name',n_pcs=30)
## 降维+聚类
sc.tl.umap(adata_con)
sc.tl.leiden(adata_con,resolution=0.5)


adata_con.write('output_path.h5ad')