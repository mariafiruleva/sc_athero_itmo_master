suppressMessages(library(argparse))
suppressMessages(library(dplyr))
library(dyno)
library(tradeSeq)
library(Seurat)
library(ggplot2)


parser <-
  ArgumentParser(description = 'Analyze diff expression across time using trade-seq')
parser$add_argument('--rda',
                    type = "character",
                    help = 'Assay name: RNA / SCT / integrated')
parser$add_argument('--traj',
                    type = "character",
                    help = 'Assay name: RNA / SCT / integrated')
parser$add_argument('--sds',
                    type = "character",
                    help = 'Path to the Seurat object in RData format')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'path to output dir')
parser$add_argument('--thr',
                    type = "integer",
                    help = 'Number of threads')
parser$add_argument('--nknots',
                    type = "integer",
                    help = 'Number of knots')

args <- parser$parse_args()
print(args)

out_dir <- args$out_dir
rda <- args$rda
thr <- args$thr
sds <- args$sds
nknots <- args$nknots
traj <- args$traj

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- thr

set.seed(1)

load(traj)
load(sds)

pseudotime <- model$sling_out$pseudotime %>%
  tibble::rownames_to_column('cell_id') %>% 
  filter(!is.na(curve1)) %>% filter(!is.na(curve2)) %>% 
  tibble::column_to_rownames('cell_id')
cell_weights <- model$sling_out$cell_weights %>%
  tibble::rownames_to_column('cell_id') %>%
  filter(curve1 > 0) %>% filter(curve2 > 0) %>% 
  tibble::column_to_rownames('cell_id')


load(rda)
counts <- as.matrix(obj@assays$RNA@counts[rownames(obj@assays$integrated), rownames(pseudotime)])


mat <- obj$GSE %>% as.data.frame() %>%
  magrittr::set_colnames('study') %>%
  tibble::rownames_to_column('cell') %>%
  filter(cell %in% colnames(counts))

rm(rda)

sce <- fitGAM(counts = counts, pseudotime=pseudotime, cellWeights = cell_weights,
              nknots = nknots, verbose = T, parallel = T, BPPARAM = BPPARAM, U=model.matrix(~ study, mat))

get_plot <- function(sce, counts, gene, test_name, out_dir) {
  plot_dir <- sprintf('%s/genes/%s', out_dir, test_name)
  plotSmoothers(sce, counts, gene)+
    ggtitle(sprintf('%s: %s', test_name, gene))+
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = T)
  }
  ggsave(sprintf('%s/%s.png', plot_dir, gene))
}

get_tbl <- function(df, test_name, out_dir) {
  tbl_dir <- sprintf('%s/tables', out_dir)
  if (!dir.exists(tbl_dir)) {
    dir.create(tbl_dir)
  }
  write.table(df, row.names = F, quote = F, file = sprintf('%s/%s.tsv', tbl_dir, test_name), sep='\t')
}


assoRes <- associationTest(sce)
startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)[1:20]

endRes <- diffEndTest(sce)
oEndRes <- order(endRes$waldStat, decreasing = TRUE)[1:20]

patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)[1:20]

## save imgs

save.image(file=sprintf('%s/trade_seq.RData', out_dir))

## save tables

get_tbl(startRes, 'start_vs_end_test', out_dir)
get_tbl(endRes, 'diff_end_test', out_dir)
get_tbl(patternRes, 'pattern_test', out_dir)


## get plots for top-20 genes

sapply(oStart, function(x) get_plot(sce, counts, gene = names(sce)[x], test_name = 'startVsEndTest', out_dir = out_dir))
sapply(oEndRes, function(x) get_plot(sce, counts, gene = names(sce)[x], test_name = 'diffEndTest', out_dir = out_dir))
sapply(oPat, function(x) get_plot(sce, counts, gene = names(sce)[x], test_name = 'patternTest', out_dir = out_dir))

