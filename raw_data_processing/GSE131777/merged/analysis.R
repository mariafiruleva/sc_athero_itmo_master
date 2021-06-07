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
  data[['percent.mito_log10']] <- log10(data[['percent.mito']])
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']])
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']])
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']])
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']])
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

filter_mito <- function(dataset, path){
  mt_dist <- as.data.frame(dataset[['scaled_mito']][[1]])
  mt_pers <- as.data.frame(dataset[['percent.mito']][[1]])
  scaled_mito_percentage <- scale(dataset[['percent.mito']][[1]])
  colnames(mt_dist) <- 'scaled_mito'
  colnames(mt_pers) <- 'percent.mito'
  filtration_coord <- get_conf_interval(dataset, 'scaled_mito')[2] * 
    attr(scaled_mito_percentage, 'scaled:scale') + 
    attr(scaled_mito_percentage, 'scaled:center')
  ggplot(mt_pers, aes(x = percent.mito)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    geom_vline(xintercept=filtration_coord, colour = "red") +
    theme(aspect.ratio = 1) +
    ggtitle('percent.mito distribution before filtration')
  ggsave(paste0(paste0(path, unique(dataset$sample)), '_mt_content_hist_before_filtration.pdf'))
  expr <- FetchData(object = dataset, vars = 'scaled_mito')
  dataset <- dataset[, which(x = expr < get_conf_interval(dataset, 'scaled_mito')[2])]
  dataset
}


## GATHERING DATA TOGETHER
    

options(future.globals.maxSize = 10000 * 1024^2)

whole <- c()
names <- c()

data <- get(load("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824239/counts.RData"))
data <- subset(x = data, features = (rowSums(as.matrix(GetAssayData(object = data, slot = "counts"))) > round(ncol(data) * 0.001, 0)))
data <- add_metadata(data)
data$sample <- "SRS4824239"
assign(paste0('data_', "SRS4824239"), data)
whole <- c(get(paste0('data_', "SRS4824239")), whole)
names <- c("SRS4824239", names) # extract SRS from /path/SRA*_SRS*.sparse.RData
data <- get(load("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824240/counts.RData"))
data <- subset(x = data, features = (rowSums(as.matrix(GetAssayData(object = data, slot = "counts"))) > round(ncol(data) * 0.001, 0)))
data <- add_metadata(data)
data$sample <- "SRS4824240"
assign(paste0('data_', "SRS4824240"), data)
whole <- c(get(paste0('data_', "SRS4824240")), whole)
names <- c("SRS4824240", names) # extract SRS from /path/SRA*_SRS*.sparse.RData
data <- get(load("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824241/counts.RData"))
data <- subset(x = data, features = (rowSums(as.matrix(GetAssayData(object = data, slot = "counts"))) > round(ncol(data) * 0.001, 0)))
data <- add_metadata(data)
data$sample <- "SRS4824241"
assign(paste0('data_', "SRS4824241"), data)
whole <- c(get(paste0('data_', "SRS4824241")), whole)
names <- c("SRS4824241", names) # extract SRS from /path/SRA*_SRS*.sparse.RData
data <- get(load("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824242/counts.RData"))
data <- subset(x = data, features = (rowSums(as.matrix(GetAssayData(object = data, slot = "counts"))) > round(ncol(data) * 0.001, 0)))
data <- add_metadata(data)
data$sample <- "SRS4824242"
assign(paste0('data_', "SRS4824242"), data)
whole <- c(get(paste0('data_', "SRS4824242")), whole)
names <- c("SRS4824242", names) # extract SRS from /path/SRA*_SRS*.sparse.RData
whole <- as.list(whole)
names(whole) <- names
rm(list=ls(pattern='data'))

if (length(whole) > 20) {
  ref <- 10
} else {
  ref <- NULL
}

## CREATE DIRECTORY FOR PLOTS
path <- paste0('plots_', "GSE131777")
dir.create(path)
path <- paste0(path, '/')

## Number of cells before

cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])


## FILTER MT CONTENT

whole <- sapply(whole, function(x) filter_mito(x, path))


## NORMALIZATION
whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  verbose = T,
  conserve.memory = T
))
whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)
whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features, 
                            verbose = FALSE)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT", 
                                        anchor.features = whole.features, verbose = FALSE, reference = ref)
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT", 
                                  verbose = FALSE)

gc()

## PCA

whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)


## TSNE

whole.integrated <- RunTSNE(whole.integrated, dims = 1:20, tsne.method = "FIt-SNE",
                            fast_tsne_path = "/nfs/home/kzajcev/FIt-SNE/bin/fast_tsne", nthreads = 4, max_iter = 2000)


## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)


## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
whole.integrated <- FindClusters(object = whole.integrated, resolution = 0.6)


## VISUALIZATION

draw_plots(path, whole.integrated)


## AVERAGING

cluster.averages <- AverageExpression(object = whole.integrated, assays = 'SCT', slot = 'data')
sapply(names(cluster.averages), 
       function(x) write.table(cluster.averages[[x]], file=paste0(x, "_clusters.tsv")))



## FINDING ANS SAVING MARKERS

whole.markers <- FindAllMarkers(object = whole.integrated,
                                assay='SCT',
                                only.pos = TRUE,
                                min.pct = 0.10,
                                test.use = 'MAST')


write.table(whole.markers, "markers.tsv", sep="\t", quote=F, row.names=F)


whole.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top50_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top100_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
top200_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)

top50_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = p_val_adj)
top100_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = p_val_adj)
top200_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = p_val_adj)

## SAVING


file_out <- paste0("GSE131777", '.RData')
save(list = c('whole.integrated', 'whole.markers', 'whole.features', 'whole.anchors'), file = file_out)

write.table(top50_log_fc, "top50_log_fc.tsv", sep="\t", quote=F, row.names=F)
write.table(top100_log_fc, "top100_log_fc.tsv", sep="\t", quote=F, row.names=F)
write.table(top200_log_fc, "top200_log_fc.tsv", sep="\t", quote=F, row.names=F)
write.table(top50_adj_pval, "top50_adj_pval.tsv", sep="\t", quote=F, row.names=F)
write.table(top100_adj_pval, "top100_adj_pval.tsv", sep="\t", quote=F, row.names=F)
write.table(top200_adj_pval, "top200_adj_pval.tsv", sep="\t", quote=F, row.names=F)
  

library(reticulate)
use_condaenv('anaconda3', required = T)
loompy <- reticulate::import('loompy')

sample <- fread(paste("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777", paste0(list.files(path = "/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777")[grepl("SRS", list.files(path = "/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777"))][1], '/sample_description.csv'), sep='/'))
adata <- sceasy:::seurat2anndata(whole.integrated, outFile=paste0(unique(sample$GSE), ".h5ad"),
                                 assay="SCT", main_layer="counts")


if (nchar(unique(sample$GSE), type = "chars", allowNA = FALSE) == 9) {
  brief <- fread(paste0('/scratch/mfiruleva/autumn/find_scRNAseqs_brief/series/',
                        paste0(paste0(stringr::str_extract(unique(sample$GSE), 'GSE[0-9]{3}'), 'nnn/'),
                               paste0(unique(sample$GSE), '/matrix/brief'))),
                 header = F, sep= '=') %>% magrittr::set_colnames(c("annot", "data"))
} else {
  brief <- fread(paste0('/scratch/mfiruleva/autumn/find_scRNAseqs_brief/series/',
                        paste0(paste0(stringr::str_extract(unique(sample$GSE), 'GSE[0-9]{2}'), 'nnn/'),
                               paste0(unique(sample$GSE), '/matrix/brief'))),
                 header = F, sep= '=') %>% magrittr::set_colnames(c("annot", "data"))
}


techs <- c()
for (file in list.files(path = "/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777")[grepl("SRS", list.files(path = "/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777"))]) {
  techs <- c(unique(fread(paste("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777", paste(file, 'sample_description.csv', sep='/'), sep = '/'))$technology), techs)
}
techs <- unique(techs)

if (techs == '10x') {
  techs <- c()
  for (file in list.files(path = "/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777")[grepl("SRS", list.files(path = "/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777"))]) {
    techs <- c(stringr::str_extract(readr::read_file(paste("/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777", paste(file, 'kallisto.sh', sep='/'), sep = '/')), "10xv[1-3]"), techs)
  }
}
techs <- unique(techs)

adata$uns["token"] = unique(sample$GSE)
adata$uns["species"] = unique(sample$scientific_name)
adata$uns["public"] = TRUE
adata$uns["curated"] = FALSE
adata$uns["sample_number"] = length(unique(whole.integrated$sample))
adata$uns["samples"] = paste(unique(whole.integrated$sample), collapse = ',')
adata$uns["study_accession"] = unique(sample$study_accession)
adata$uns["title"] = filter(brief,  annot=='!Series_title')$data
adata$uns["description"] = filter(brief,  annot=='!Series_summary')$data
adata$uns["design"] = filter(brief,  annot=='!Series_overall_design')$data
adata$uns["geo"] = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", unique(sample$GSE))
adata$uns["pubmed"] = paste0("https://www.ncbi.nlm.nih.gov/pubmed/", filter(brief,  annot=='!Series_pubmed_id')$data)
adata$uns['status'] = filter(brief,  annot=='!Series_status')$data
adata$uns["technology"] = paste(techs,collapse=",")
adata$uns["cells"] = ncol(whole.integrated)
adata$uns["expType"] = "counts"

adata$write(paste0(unique(sample$GSE), ".h5ad"), compression = "gzip")

## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)



