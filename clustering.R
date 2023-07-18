library(Seurat)
library(harmony)
library(dplyr)
library(SeuratObject)


preprocessing = function(data,filter_cells, n_genes_by_counts_max = 60000, n_genes_by_counts_min = 200, pct_counts_mt_min = 10){
    data[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
    data[["percent.rb"]] = PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
    data[["percent.hb"]] = PercentageFeatureSet(pbmc, pattern = "^HB[APS]")
    data[["percent.hsp"]] = PercentageFeatureSet(pbmc, pattern = "^HSP")
    data = subset(data, subset = nFeature_RNA > n_genes_by_counts_min & 
                  nFeature_RNA < n_genes_by_counts_max & 
                  percent.mt < pct_counts_mt_min)
    data = data[,setdiff(colnames(data), filter_cells)]
    return(data)
}


## 单个样本预处理
datalist = lapply(X = datalist, FUN = function(x) {
    x = preprocessing(x, c())
})

## 多个样本合起来
gene_names = rownames(datalist[[0]])
for(data in datalist[2:length(datalist)]){
    gene_names = intersect(gene_names, rownames(data))
}

combined_count = GetAssayData(datalist[[0]], slot = "count", assay = 'RNA')
combined_count = combined_count[gene_names, ]
group = rep('g1',ncol(combined_count))
i = 2
for(data in datalist[2:length(datalist)]){
    current_count = GetAssayData(data, slot = "count", assay = 'RNA')
    combined_count = cbind(combined_count, current_count)
    group = c(group, rep(paste0('g',i),ncol(current_count)))
}

## 多个样本合起来标化
combined = CreateSeuratObject(counts = combined_count) %>%
       Seurat::NormalizeData(verbose = FALSE) %>%
       FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
       ScaleData(verbose = FALSE,vars.to.regress = c('nFeature_RNA','nCount_RNA','percent.mt') %>% 
       RunPCA(pc.genes = combined@var.genes, npcs = 30, verbose = FALSE)
combined$batch = combined
## 多个样本去批次效应
combined = RunHarmony(combined, "batch", plot_convergence = FALSE))
## 降维+聚类
combined = RunUMAP(combined, reduction='harmony', dims = 1:30) %>% 
        FindNeighbors(dims=1:30,k.param = 25,reduction = 'harmony') %>% 
        FindClusters(resolution=0.5,reduction = 'harmony')


saveRDS(combined, 'output_path.rds')                 
                 