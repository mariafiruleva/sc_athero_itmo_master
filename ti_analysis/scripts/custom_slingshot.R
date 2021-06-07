suppressMessages(library(argparse))
suppressMessages(library(babelwhale))
suppressMessages(library(dyno))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))

parser <-
  ArgumentParser(description = 'Trajectory analysis using dyno wrapper')
parser$add_argument('--assay_data',
                    type = "character",
                    help = 'Assay name: RNA / SCT / integrated')
parser$add_argument('--assay_counts',
                    type = "character",
                    help = 'Assay name: RNA / SCT / integrated')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to the Seurat object in RData format')
parser$add_argument('--ident',
                    type = "character",
                    help = 'Ident name, e.g., integrated_snn_res.0.6')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'path to output dir')
parser$add_argument('--clusters',
                    type = "integer", nargs='+',
                    help = 'Clusters ids')
parser$add_argument('--start_clus',
                    type = "integer",
                    help = 'ID of start cluster')

## SET VARIABLES

args <- parser$parse_args()
print(args)

## DEFINE SOME FUCNTIONS

get_data <- function(path, ident, clusters, start_clus, assay_counts='SCT', assay_data='integrated') {
  obj <- get(load(path))
  Idents(obj) <- ident
  obj <- subset(obj, idents = clusters)
  wrap_data <- wrap_expression(
    counts = t(as.matrix(obj@assays[[assay_counts]]@counts)),
    expression = t(as.matrix(obj@assays[[assay_data]]@data))
  )
  wrap_data <- add_prior_information(
    wrap_data,
    groups_id=cbind(colnames(obj),obj[[ident]][[1]]) %>%
      as.data.frame() %>%
      magrittr::set_colnames(c('cell_id', 'group_id')),
    start_id = sample(colnames(subset(obj, idents = start_clus)), 1)
  )
  wrap_data <- add_dimred(
    wrap_data,
    dimred = as.matrix(obj@reductions$umap@cell.embeddings)
  )
  wrap_data <- add_grouping(
    wrap_data,
    obj[[ident]][[1]]
  )
  wrap_data
}

get_branching_point <- function(model, obj) {
  branching_milestone <- model$milestone_network %>%
    group_by(from) %>% filter(n() > 1) %>% pull(from) %>% dplyr::first()
  branch_point <- calculate_branching_point_feature_importance(model,
                                                               expression_source=obj$expression,
                                                               milestones_oi = branching_milestone)
  branch_point
}

get_hmap <- function(x, branch_feature_importance, model, wrap_data) {
  plot_heatmap(model, expression_source = wrap_data$expression, features_oi =
                 branch_feature_importance %>%
                 filter(from == x[1]) %>%
                 filter(to == x[2]) %>%
                 top_n(50, importance) %>%
                 pull(feature_id) %>%
                 unique() %>% as.character())+
    ggtitle(sprintf('top-50, from %s to %s', x[1], x[2]))
  ggsave(sprintf('plots/hmaps/top50_from_%s_to-%s.png', x[1], x[2]))
  plot_heatmap(model, expression_source = wrap_data$expression, features_oi =
                 branch_feature_importance %>%
                 filter(from == x[1]) %>%
                 filter(to == x[2]) %>%
                 top_n(100, importance) %>%
                 pull(feature_id) %>%
                 unique() %>% as.character())+
    theme(text = element_text(size=12))+
    ggtitle(sprintf('top-100, from %s to %s', x[1], x[2]))
  ggsave(sprintf('plots/hmaps/top100_from_%s_to-%s.png', x[1], x[2]),
         width=18, height=18)
}

save_tables <- function(model, branch_feature_importance, branch_point) {
  write.table(model$progressions,
              file='tables/progression.csv', quote = F, row.names = F, sep=',')
  write.table(model$milestone_network,
              file='tables/milestone_network.csv', quote = F, row.names = F, sep=',')
  write.table(branch_feature_importance,
              file='tables/branch_feature_importance.csv', quote = F, row.names = F, sep=',')
  write.table(branch_point,
              file='tables/branch_point.csv', quote = F, row.names = F, sep=',')
}

## PREPARE THE OBJECT

wrap_data <- get_data(args$data, args$ident, args$clusters, args$start_clus)

## INFER TRAJECTORY

set.seed(1)
setwd(args$out_dir)

parameters <- list(cluster_method = "pam", ndim = 20L, shrink = 1L, reweight = TRUE,
                   reassign = TRUE, thresh = 0.001, maxit = 10L, stretch = 2L,
                   smoother = "smooth.spline", shrink.method = "cosine")
priors <- list()
priors$start_id <- wrap_data$prior_information$start_id

model <- run_fun(
  expression = wrap_data$expression,
  priors = priors,
  parameters = parameters,
  verbose = T
)

## UMAP

dimred <- dyndimred::dimred_umap(wrap_data$expression)

## CALCULATE SOME FEATURES: FEATURE IMPORTANCE, BRANCHING POINT

branch_feature_importance <- calculate_branch_feature_importance(model, expression_source=wrap_data$expression)
branch_point <- get_branching_point(model, wrap_data)

## VISUALIZATION

plot_dimred(
  model,
  dimred = dimred,
  expression_source = wrap_data$expression,
  color_density = "grouping",
  alpha_cells = 0.7,
  grouping = wrap_data$prior_information$groups_id,
  label_milestones = T,
  hex_cells = F
)+
  theme(aspect.ratio = 1, legend.position = "none")+
  ggtitle(sprintf('clusters, %s', args$ident))
ggsave('plots/trajectory_clusters.png')

dir.create('plots/hmaps')
apply(model$milestone_network %>% select(from, to), 1, function(x) get_hmap(x, branch_feature_importance, model, wrap_data))

## SAVE OUTPUT: RDATA

save(list = c('args', 'wrap_data', 'model', 'dimred', 'branch_feature_importance', 'branch_point'),
     file = 'rdata/trajectory.RData')

## SAVE OUTPUT: TABLES

# save_tables(model, branch_feature_importance, branch_point)