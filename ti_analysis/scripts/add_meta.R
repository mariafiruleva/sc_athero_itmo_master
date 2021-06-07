suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(functools))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(Matrix))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))


set.seed(1)

parser <-
  ArgumentParser(description = 'Reanalyze the prepared seurat rdata object')
parser$add_argument('--in_rda',
                    type = "character",
                    help = 'path to the rdata object')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'object name')
parser$add_argument('--traj',
                    type = "character",
                    help = 'path to rda with trajectory results')
args <- parser$parse_args()
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
  seurat_object
}


## ADD METADATA USING TRAJECTORY ANALYSIS RESULTS

obj <- add_meta(args$traj, obj)

## SAVE: RDATA

setwd(args$out_dir)

save(list = c('obj'), file = "object_with_meta.RData")

## SAVE: h5

library(reticulate)
use_condaenv('anaconda3', required = T)
sceasy:::seurat2anndata(obj, outFile="object_with_meta.h5ad", assay="SCT", main_layer="counts")