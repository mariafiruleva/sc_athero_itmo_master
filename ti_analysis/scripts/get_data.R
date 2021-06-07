suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(functools))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(Matrix))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))


set.seed(1)

args <- list()
args$in_rda <- 'target/data/target.RData'
args$out_dir <- 'out/slingshot/rdata'
args$traj <- 'out/slingshot/rdata/trajectory.RData'

print(args)

load(args$in_rda)

add_meta <- function(trajectory, seurat_object) {
  load(trajectory)
  dimred <- as.data.frame(dimred) %>%
    tibble::rownames_to_column(var='cell_id') %>%
    arrange(match(cell_id, colnames(seurat_object)))
  model$milestone_percentages <- model$milestone_percentages %>% arrange(cell_id, -percentage) %>%
    filter(!duplicated(cell_id)) %>%
    arrange(match(cell_id, colnames(seurat_object)))
  old_meta <- cbind(names(wrap_data$grouping), wrap_data$grouping) %>%
    magrittr::set_colnames(c('cell_id', 'cluster')) %>%
    as.data.frame() %>% 
    arrange(match(cell_id, colnames(seurat_object)))
  seurat_object$DYNO_UMAP_1 <- dimred$comp_1
  seurat_object$DYNO_UMAP_2 <- dimred$comp_2
  seurat_object$lineage <- model$milestone_percentages$milestone_id
  seurat_object$old_clusters <- old_meta$cluster
  seurat_object
}

## PCA

obj <- RunPCA(obj, verbose = FALSE)

## TSNE

whole.integrated <- RunTSNE(obj, dims = 1:20, tsne.method = "FIt-SNE",
                            fast_tsne_path = "/opt/conda/pkgs/fit-sne-1.1.0-h3ddc34e_0/bin/fast_tsne", nthreads = 4, max_iter = 2000)

## UMAP

obj <- RunUMAP(obj, dims = 1:20)


## CLUSTERING

obj <- FindNeighbors(object = obj, dims = 1:20)
obj <- FindClusters(object = obj, resolution = seq(0.2, 1, 0.2))

## ADD METADATA USING TRAJECTORY ANALYSIS RESULTS

obj <- add_meta(args$traj, obj)

## SAVE: RDATA

setwd(args$out_dir)

save(list = c('obj'), file = "object.RData")

## VISUALIZATION

analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  cluster.averages <- AverageExpression(object = object, assays = 'SCT', slot = 'data')
  dir.create(sprintf("markers/%s", ident), recursive=T)
  sapply(names(cluster.averages),
         function(x) write.table(cluster.averages[[x]], file=sprintf('markers/%s/%s_clusters.tsv', ident, x)))

  ## FINDING ANS SAVING MARKERS

  whole.markers <- FindAllMarkers(object = object,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'MAST')
  write.table(whole.markers, sprintf('markers/%s/markers.tsv', ident), sep="\t", quote=F, row.names=F)
}


## SAVING: h5

library(reticulate)
sceasy:::seurat2anndata(obj, outFile="object.h5ad", assay="SCT", main_layer="counts")

## MARKERS

idents <- c('integrated_snn_res.0.2', 'integrated_snn_res.0.4', 'integrated_snn_res.0.6', 'integrated_snn_res.0.8', 'integrated_snn_res.1')
dir.create('markers')
sapply(idents, function(ident) analyze_object(object = obj, ident = ident))
