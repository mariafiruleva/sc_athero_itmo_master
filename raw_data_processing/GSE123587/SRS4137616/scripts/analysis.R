suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(functools))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(reticulate))

setwd("/mnt/tank/scratch/mfiruleva/scn/data/GSE123587/SRS4137616")
set.seed(1)

## FUNCTIONS


add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mito_log10']] <- log10(data[['percent.mito']] + 1)
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']] + 1)
  data[['scaled_mito']] <- scale(percent.mito)
  data[['scaled_nCount_RNA']] <- scale(data[['nCount_RNA_log10']])
  attr(data$scaled_nCount_RNA, "scaled:center") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:scale") <- NULL
  attr(data$scaled_mito, "scaled:center") <- NULL
  attr(data$scaled_mito, "scaled:scale") <- NULL
  data
}

draw_plots <- function(path, data) {
  VlnPlot(
    data,
    features = "nFeature_RNA",
    pt.size = 0.1
    ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_features.pdf'))
  VlnPlot(
    data,
    features = "nCount_RNA",
    pt.size = 0.1
  ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_umi.pdf'))
  VlnPlot(
    data,
    features = "percent.mito",
    pt.size = 0.1
  ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_mt.pdf'))
  
  
  
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mito") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_mt_plot.pdf'))
  
  FeatureScatter(data, feature1 = "nCount_RNA_log10", feature2 = "nFeature_RNA_log10") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_features_log10_plot.pdf'))
  
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_features_plot.pdf'))
  
  ElbowPlot(data, ndims = 50) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'elbow_plot.pdf'))
}

get_conf_interval <- function(dataset, parameter) {
  left <- mean(dataset[[parameter]][[1]]) - qnorm(0.975)
  right <- mean(dataset[[parameter]][[1]]) + qnorm(0.975)
  return(c(left, right))
}


## GATHERING DATA TOGETHER
    
options(future.globals.maxSize = 8000 * 1024^2)

whole <- get(load("/mnt/tank/scratch/mfiruleva/scn/data/GSE123587/SRS4137616/counts.RData"))
whole <- subset(x = whole, features = (rowSums(as.matrix(GetAssayData(object = whole, slot = "counts"))) > round(ncol(whole) * 0.001, 0)))
whole <- add_metadata(whole)




## CREATE DIRECTORY FOR PLOTS
path <- "plots_SRS4137616/"
dir.create(path)

## Number of cells before

cells.before <- dim(GetAssayData(object = whole, slot = "counts"))[2]


## FILTER MT CONTENT

mt_dist <- as.data.frame(whole[['scaled_mito']][[1]])
colnames(mt_dist) <- 'scaled_mito'
ggplot(mt_dist, aes(x = scaled_mito)) +
  geom_histogram(color = "black", fill = "white", bins = 30) +
  geom_vline(xintercept=get_conf_interval(whole, 'scaled_mito')[2], colour = "red") +
  ggtitle('scaled_mito distribution before filtration')
ggsave(paste0(path, 'mt_content_hist_before_filtration.pdf'))
whole <-
  subset(
    x = whole,
    subset = scaled_mito < get_conf_interval(whole, 'scaled_mito')[2]
  )


## NORMALIZATION

whole <-
  SCTransform(
    whole,
    ncells=min(100000, ncol(whole)),
    vars.to.regress = c("percent.mito"),
    verbose = T,
    conserve.memory = T
  )





## PCA
gc()

whole <- RunPCA(object = whole, features = VariableFeatures(object = whole), npcs=50)


## TSNE

whole <-
  RunTSNE(whole, dims = 1:20, tsne.method = "FIt-SNE",
          fast_tsne_path = "/opt/conda/pkgs/fit-sne-1.1.0-h3ddc34e_0/bin/fast_tsne", nthreads = 4, max_iter = 2000)


## UMAP

whole <- RunUMAP(whole, dims = 1:20)


## CLUSTERING

whole <- FindNeighbors(object = whole, dims = 1:20)
whole <- FindClusters(object = whole, resolution = c(0.2, 0.4, 0.6, 0.8, 1))


## VISUALIZATION

draw_plots(path, whole)


## AVERAGING

cluster.averages <- AverageExpression(object = whole, assays = 'SCT', slot = 'data')
sapply(names(cluster.averages), 
       function(x) write.table(cluster.averages[[x]], file=paste0(x, "_clusters.tsv")))


## FINDING ANS SAVING MARKERS

analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  out_dir <- paste0('markers/', ident)
  dir.create(out_dir, recursive = T)
  whole.markers <- FindAllMarkers(object = object,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'MAST')
  write.table(whole.markers, paste(out_dir, "markers.tsv", sep = '/'), sep="\t", quote=F, row.names=F)
  top50_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
  top100_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
  top200_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)

  top50_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = p_val_adj)
  top100_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = p_val_adj)
  top200_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = p_val_adj)

  write.table(top50_log_fc, paste0(out_dir, "/top50_log_fc.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top100_log_fc, paste0(out_dir, "/top100_log_fc.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top200_log_fc, paste0(out_dir, "/top200_log_fc.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top50_adj_pval, paste0(out_dir, "/top50_adj_pval.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top100_adj_pval, paste0(out_dir, "/top100_adj_pval.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top200_adj_pval, paste0(out_dir, "/top200_adj_pval.tsv"), sep="\t", quote=F, row.names=F)
}


idents <- paste0('SCT_snn_res.', c(0.2, 0.4, 0.6, 0.8, 1))
sapply(idents, function(ident) analyze_object(object = whole, ident = ident))




## SAVING


save('whole', file = "SRS4137616.RData")
  

library(reticulate)

sceasy:::seurat2anndata(whole, outFile="SRS4137616.h5ad",
                        assay="SCT", main_layer="counts")

## Number of cells after

cells.after <- length(colnames(x = whole))
print(paste0("cells.before:",cells.before))
print(paste0("cells.after:",cells.after))
print(paste0("cell.diff:", cells.before-cells.after))




