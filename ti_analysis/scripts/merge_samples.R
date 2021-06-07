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

set.seed(1)


## GATHERING DATA TOGETHER


options(future.globals.maxSize = 4e4 * 1024^2)

get_whole_obj <- function(pathes) {
  objects <- lapply(pathes, function(x) SplitObject(get(load(x)), split.by = 'GSE'))
  objects <- do.call(c, unlist(objects, recursive=F))
  names(objects) <- sapply(objects, function(x) unique(x$GSE))
  objects
}

files <- list.files('target', full.names = T)

print(files)

whole <- get_whole_obj(files)

## NORMALIZATION

whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  verbose = F,
  conserve.memory = T
))

## INTEGRATION

whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 3000, assay=rep('SCT', length(names(whole))))
whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features, verbose = F)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT", anchor.features = whole.features, verbose = F)
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
                                  verbose = FALSE)

## PCA
gc()

whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)


## TSNE

whole.integrated <- RunTSNE(whole.integrated, dims = 1:30, tsne.method = "FIt-SNE",
                            fast_tsne_path = "/opt/conda/pkgs/fit-sne-1.1.0-h3ddc34e_0/bin/fast_tsne", nthreads = 4, max_iter = 2000)


## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:30)


## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:30)
whole.integrated <- FindClusters(object = whole.integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1))


## SAVING


library(reticulate)

sceasy:::seurat2anndata(whole.integrated, outFile="athero_merged.h5ad",
                        assay="SCT", main_layer="counts")

save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = "athero_merged.RData")

## AVERAGING

cluster.averages <- AverageExpression(object = whole.integrated, assays = 'SCT', slot = 'data')
sapply(names(cluster.averages), 
       function(x) write.table(cluster.averages[[x]], file=paste0(x, "_clusters.tsv")))



## FINDING ANS SAVING MARKERS

analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  if (length(levels(object)) == 1) {
    return(message(sprintf('%s: since only one cluster was identified, markers can not be found', ident)))
  }
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

idents <- paste0('integrated_snn_res.', c(0.2, 0.4, 0.6, 0.8, 1))
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))