library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Matrix)
library(dplyr)
library(grid)
library(gridExtra)
library(pheatmap)
library(gtable)
library(plyr)
library(ggsignif)
library(MASS)
library(patchwork)
library(ComplexHeatmap)

load('/mnt/tank/scratch/mfiruleva/scn/private_data/bat_stoyan/stoyan.RData')

getPalette.1 <- colorRampPalette(brewer.pal(9, "Set1"))

## X_A

str(whole.integrated$integrated_snn_res.1)

whole.integrated$custom_clusters <- whole.integrated$integrated_snn_res.1
## num lbls

whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '0', '0', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '11', '1', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '2', '2', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '9', '3', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '1', '4', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '8', '5', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '7', '6', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '12', '7', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '13', '8', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '10', '9', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '5' | whole.integrated$integrated_snn_res.1 == '6' | whole.integrated$integrated_snn_res.1 == '14', '10', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '4' | whole.integrated$integrated_snn_res.1 == '16', '11', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '3' | whole.integrated$integrated_snn_res.1 == '15', '12', whole.integrated$custom_clusters)


sorted_labels <- 0:12

whole.integrated$custom_clusters <- factor(x = whole.integrated$custom_clusters, levels = sorted_labels)
Idents(whole.integrated) <- 'custom_clusters'
new_labels <- c('0 Neutrophils', expression(~1~Myeloid~cells~italic(Zeb2)^{"hi"}), expression(~2~Monocytes~italic(Ly6C)^{"low"}),
                                            expression(~3~Monocytes~italic(Ly6C)^{"int"}), expression(~4~Monocytes~italic(Ly6C)^{"hi"}),
                                            '5 Matrix Macrophages', '6 Macrophages M2-like', expression(7~Macrophages~italic(Lpl)^{"hi"}), expression(8~Macrophages~italic(Plin2)^{"hi"}),
                '9 Dendritic cells', '10 T cells', '11 B cells', '12 NK cells')

whole.integrated$genotype <- factor(x = whole.integrated$genotype, levels = c("WT", "KO"))
plt <- DimPlot(whole.integrated, split.by='genotype', pt.size=0.25)+
  scale_color_manual(values=getPalette.1(length(unique(whole.integrated$custom_clusters))), labels = new_labels)+
  theme_bw(base_size=11)+
  theme(legend.text.align = 0, legend.key.size=unit(0.2, "in"), aspect.ratio = 1, plot.margin=grid::unit(c(0,0,0.2,0), "in"),
        #legend.position = 'none',
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_legend(ncol = 4, override.aes = list(color = c('black'))))+
  scale_fill_continuous(guide="legend",breaks=seq(0.2,0.8,by=0.1))

LabelClusters(plt, id = "ident", size=5, fontface = 'bold', repel=F)

ggsave("clustering_total.png", width = 8, height = 4, dpi=600, units='in')

ggsave("clustering_total.svg", width = 11, height = 5, dpi=600, units='in', device = 'svg')

## X_B (add asterics)

df <- cbind(as.character(whole.integrated$genotype), as.character(whole.integrated$custom_clusters)) %>% as.data.frame() %>% 
  magrittr::set_colnames(c('genotype', 'cluster'))

df$cluster <- factor(x = df$cluster, levels = as.character(sort(as.numeric(levels(df$cluster)))))
df$genotype <- factor(x = df$genotype, levels = c("WT", "KO"))

p <- df %>% group_by(cluster, genotype) %>% dplyr::summarise(n = n())

p <- p %>%
  group_by(genotype) %>%
  mutate(total = sum(n)) %>%
  group_by(cluster, add=TRUE) %>%
  mutate(per=round(100*n/total,2))

ggplot(p,aes(x=cluster,y=per, fill=genotype))+ 
  geom_bar(stat = 'identity', position=position_dodge(width = 0.5), width = 0.45) + 
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0, legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  geom_segment(aes(x = 2, y = 8, xend = 2, yend = 6),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 6, y = 13, xend = 6, yend = 11),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 7, y = 13, xend = 7, yend = 11),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 8, y = 10, xend = 8, yend = 8),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 9, y = 10, xend = 9, yend = 8),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  ylab('cells per clutser, %')


ggsave("hist.png", width = 4, height = 4, dpi=600, units='in')
ggsave("hist.svg", width = 4, height = 4, dpi=600, units='in')


## X_F


get_expr <- function(object, gene_set, slot='data', assay='SCT') {
  av <- numeric(ncol(object))
  zz <- which(tolower(rownames(GetAssayData(object, slot = slot, assay = assay))) %in% tolower(gene_set))
  object@assays$SCT@data[zz, ]
}

source('split_violin.R')

plot_vln <- function(object, gene_set, pw_name, reduction="umap", assay='SCT', slot='data') {
  red <- cbind(get_expr(object, gene_set)) %>% as.data.frame() %>% 
    magrittr::set_colnames(c('expression'))
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c("WT", "KO"))
  red$cluster <- object$custom_clusters
  target_clusters <- c(2, 4, 3, 6, 5, 7, 8)
  red <- red %>% dplyr::filter(cluster %in% target_clusters)
  ggplot(data=red, aes(x=cluster, y=expression, fill=genotype)) + theme_bw(base_size = 8) + 
    theme(panel.spacing = unit(0, "lines"),
          legend.position = 'none',
          legend.title = element_blank(),
          aspect.ratio = 0.2,
          # axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "italic", size=6, margin=margin(0,0,0,0)),
          plot.margin=grid::unit(c(0,0,0,0), "in")) + 
    geom_split_violin(scale="width", alpha=0.7) +
    scale_fill_brewer(palette='Set1', guide=guide_legend(ncol=2), direction=-1) +
    ylab(NULL) + xlab(NULL)+ggtitle(gene_set)
  ggsave(sprintf("%s_vln.png", gene_set), width = 3, height = 0.8, dpi=600, units='in')
}
plot_vln(whole.integrated, 'Plin2')
plot_vln(whole.integrated, 'Lpl')
plot_vln(whole.integrated, 'Cd36')
plot_vln(whole.integrated, 'Trem2')
plot_vln(whole.integrated, 'Fabp4')
plot_vln(whole.integrated, 'Fabp5')

## trajectory
load('/mnt/tank/scratch/mfiruleva/scn/private_data/bat_stoyan/trajectory/resolution1.0/macs/slingshot/rdata/trajectory.RData')

new_traj_labels <- c('0 Neutrophils', expression(~1~Myeloid~cells~italic(Zeb2)^{"hi"}), expression(~2~Monocytes~italic(Ly6C)^{"low"}),
                expression(~3~Monocytes~italic(Ly6C)^{"int"}), expression(~4~Monocytes~italic(Ly6C)^{"hi"}),
                '5 Matrix Macrophages', '6 Macrophages M2-like', expression(7~Macrophages~italic(Lpl)^{"hi"}), expression(8~Macrophages~italic(Plin2)^{"hi"}),
                '9 Dendritic cells', '10 T cells', '11 B cells', '12 NK cells')

new_traj_labels <- c(expression(~4~Monocytes~italic(Ly6C)^{"hi"}), expression(7~Macrophages~italic(Lpl)^{"hi"}), expression(8~Macrophages~italic(Plin2)^{"hi"}),
                     expression(~2~Monocytes~italic(Ly6C)^{"low"}), '6 Macrophages M2-like', '5 Matrix Macrophages', expression(~3~Monocytes~italic(Ly6C)^{"int"}))

obj$cluster$grouping <- as.factor(obj$cluster$grouping)
obj$cluster$grouping <- ifelse(obj$cluster$grouping == 1, '4',
                               ifelse(obj$cluster$grouping == 9, '3',
                                      ifelse(obj$cluster$grouping == 8, '5',
                                             ifelse(obj$cluster$grouping == 7, '6',
                                                    ifelse(obj$cluster$grouping == 12, '7',
                                                    ifelse(obj$cluster$grouping == 13, '8', '2'))))))
obj$cluster$grouping <- as.factor(obj$cluster$grouping)
new_traj_labels <- c(expression(~2~Monocytes~italic(Ly6C)^{"low"}), expression(~3~Monocytes~italic(Ly6C)^{"int"}),
                     expression(~4~Monocytes~italic(Ly6C)^{"hi"}), '5 Matrix Macrophages', '6 Macrophages M2-like',
                     expression(7~Macrophages~italic(Lpl)^{"hi"}), expression(8~Macrophages~italic(Plin2)^{"hi"}))

tra <- plot_dimred(
  model,
  dimred = dimred,
  expression_source = obj$cluster$expression,
  color_density = "grouping",
  alpha_cells = 0.7,
  size_cells = 1,
  size_trajectory = 1.5,
  grouping = obj$cluster$grouping,
  label_milestones = F,
)+
  scale_fill_manual(values=getPalette.1(length(unique(whole.integrated$custom_clusters)))[3:9], labels = NULL)+
  scale_color_manual(values=getPalette.1(length(unique(whole.integrated$custom_clusters)))[3:9],labels = new_traj_labels)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=3)))+
  theme_void(base_size = 11)+
  theme(aspect.ratio = 1, legend.title = element_blank(),
        legend.text.align = 0, legend.key.size=unit(0.2, "in"), plot.margin=grid::unit(c(0,0,0.2,0), "in"), 
        legend.text=element_text(size=11))


print(tra, vp=viewport(angle=-20))
tra
ggsave("trajectory.png", width = 8, height = 6, dpi=600, units='in')
ggsave("trajectory.svg", width = 8, height = 6, dpi=600, units='in', device = 'svg')


## PATHWAY

load('/scratch/mfiruleva/autumn/fgsea/curatedSymbol.rdata')

get_expr <- function(object, gene_set, slot='data', assay='SCT') {
  av <- numeric(ncol(object))
  zz <- which(tolower(rownames(GetAssayData(object, slot = slot, assay = assay))) %in% tolower(gene_set))
  geneExp <- as.matrix(log2(object@assays$SCT@data[zz, ] + 1))
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  av <- av + colSums(geneExp) / length(gene_set)
  av
}


plot_target_pw <- function(object, gene_set, pw_name, reduction="umap", assay='SCT', slot='data') {
  red <- cbind(get_expr(object, gene_set), object@reductions[['umap']]@cell.embeddings) %>% as.data.frame() %>% 
    magrittr::set_colnames(c('expression', paste0(reduction, 1), paste0(reduction, 2)))
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c('WT', 'KO'))
  ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.15)+
    theme_bw(base_size=8) +
    facet_grid(.~genotype)+
    scale_color_gradientn(colours=c("darkblue", "blue", "grey", "red", "darkred"),
                          breaks = c(0, floor(max(red$expression) * 100) / 100),
                          rescaler = function(x, from) {
      res <- numeric(length(x))
      res[x >= 1 & !is.na(x)] <- 1
      res[x < 1 & !is.na(x)] <- (x[x < 1 & !is.na(x)] + 1) / 2
      res[x < -1 & !is.na(x)] <- 0
      res[is.na(x)] <- NA
      res
    })+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          plot.margin=grid::unit(c(0,0,0,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gsub('_', ' ', pw_name))
  ggsave(sprintf("pw/%s.png", pw_name), width = 4.25, height = 2, dpi=600, units='in')
}


pws <- jsonlite::fromJSON('pathways.json')

plot_target_pw(whole.integrated, pws$KEGG_GLYCEROLIPID_METABOLISM, 'KEGG_GLYCEROLIPID_METABOLISM')
plot_target_pw(whole.integrated, pws$KEGG_GLYCEROPHOSPHOLIPID_METABOLISM, 'KEGG_GLYCEROPHOSPHOLIPID_METABOLISM')
plot_target_pw(whole.integrated, pws$PID_LYSOPHOSPHOLIPID_PATHWAY, 'PID_LYSOPHOSPHOLIPID_PATHWAY')
plot_target_pw(whole.integrated, pws$HALLMARK_FATTY_ACID_METABOLISM, 'HALLMARK_FATTY_ACID_METABOLISM')

## heatmap
set.seed(1)
library(data.table)
whole.markers <- fread('custom_markers.tsv')

top5 <- whole.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% group_by(gene) %>% 
  filter(avg_logFC == max(avg_logFC)) %>% 
  distinct()
top10 <- whole.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% group_by(gene) %>% 
  filter(avg_logFC == max(avg_logFC)) %>% 
  distinct()
DoHeatmap(whole.integrated, features = top5$gene) + NoLegend()

## add cluster

library(ComplexHeatmap)
library(dplyr)

expr <- whole.integrated@assays$integrated@scale.data[
  top10$gene[top10$gene %in% rownames(whole.integrated@assays$integrated@scale.data)], ] %>% as.matrix() # change expr values ; now integrated

cluster_anno <- whole.integrated$custom_clusters

row_anno <- as.data.frame(top10$cluster) %>% mutate(gene = top10$gene) %>%
  filter(gene %in% rownames(expr)) %>% 
  tibble::column_to_rownames(var = 'gene') %>% 
  magrittr::set_colnames(c('V1'))
quantile(expr, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1, 0, 2), c("#FF00FF", "black", "#FFFF00"))

Heatmap(expr, name = "expr", column_gap = unit(0.5, "mm"),
        cluster_columns = F, cluster_rows = F, column_split = factor(cluster_anno), row_split = factor(row_anno$V1),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = getPalette.1(length(unique(whole.integrated$custom_clusters)))))),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = getPalette.1(length(unique(whole.integrated$custom_clusters)))))),
        show_column_names = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 6)
)

## change raster? use_raster = TRUE, raster_quality = 4 -- rm



## zoom macs

get_expr <- function(object, gene_set, slot='data', assay='SCT') {
  av <- numeric(ncol(object))
  zz <- which(tolower(rownames(GetAssayData(object, slot = slot, assay = assay))) %in% tolower(gene_set))
  geneExp <- as.matrix(object@assays$SCT@data[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  av <- av + colSums(geneExp) / length(gene_set)
  av
}


plot_target_pw <- function(object, gene_set, pw_name, reduction="umap", assay='SCT', slot='data') {
  red <- cbind(get_expr(object, gene_set), object@reductions[['umap']]@cell.embeddings) %>% as.data.frame() %>% 
    magrittr::set_colnames(c('expression', paste0(reduction, 1), paste0(reduction, 2)))
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c('WT', 'KO'))
  ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.15)+
    theme_bw(base_size=8) +
    facet_grid(.~genotype)+
    scale_color_gradientn(colours=c("darkblue", "blue", "grey", "red", "darkred"),
                          breaks = c(0, floor(max(red$expression) * 100) / 100),
                          rescaler = function(x, from) {
                            res <- numeric(length(x))
                            res[x >= 1 & !is.na(x)] <- 1
                            res[x < 1 & !is.na(x)] <- (x[x < 1 & !is.na(x)] + 1) / 2
                            res[x < -1 & !is.na(x)] <- 0
                            res[is.na(x)] <- NA
                            res
                          })+
    # xlim(NA, 5)+
    # ylim(0, NA)+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          plot.margin=grid::unit(c(0,0,0,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gsub('_', ' ', pw_name))
  ggsave(sprintf("macs/pw/%s.png", pw_name), width = 4.25, height = 2, dpi=600, units='in')
}


pws <- jsonlite::fromJSON('pathways.json')

plot_target_pw(whole.integrated, pws$KEGG_GLYCEROLIPID_METABOLISM, 'KEGG_GLYCEROLIPID_METABOLISM')
plot_target_pw(whole.integrated, pws$KEGG_GLYCEROPHOSPHOLIPID_METABOLISM, 'KEGG_GLYCEROPHOSPHOLIPID_METABOLISM')
plot_target_pw(whole.integrated, pws$PID_LYSOPHOSPHOLIPID_PATHWAY, 'PID_LYSOPHOSPHOLIPID_PATHWAY')
plot_target_pw(whole.integrated, pws$HALLMARK_FATTY_ACID_METABOLISM, 'HALLMARK_FATTY_ACID_METABOLISM')

## macs : genes

plot_target_gene <- function(object, gene, reduction="umap", assay='SCT', slot='data') {
  data <- GetAssayData(object, slot = slot, assay = assay)
  red <- object@reductions[[reduction]]@cell.embeddings
  red <- as.data.frame(red)
  colnames(red) <- paste0(reduction, 1:ncol(red))
  genes_indexes <- which(rownames(data) %in% gene[[1]])
  expression_signaling <- data[genes_indexes,]
  red$expression <- expression_signaling
  #red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  #red$genotype <- factor(x = red$genotype, levels = c('WT', 'KO'))
  red$cluster <- object$integrated_snn_res.0.2
  onlyBreak <- floor(max(red$expression))
  ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.05)+
    theme_bw(base_size=11) +
    # facet_grid(.~genotype)+
    scale_color_gradientn(colours=c("grey", "red", "red3"),breaks=c(0, onlyBreak),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    #xlim(NA, 5)+
    #ylim(0, NA)+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "italic"),
          plot.margin=grid::unit(c(0,0,0.2,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gene)
  ggsave(sprintf("%s.png", gene), width = 4, height = 2.2, dpi=600, units='in')
}

plot_target_gene(whole.integrated, 'Mrc1')
plot_target_gene(whole.integrated, 'Clec10a')
plot_target_gene(whole.integrated, 'Mki67')
plot_target_gene(whole.integrated, 'Ccna2')
plot_target_gene(whole.integrated, 'Top2a')

plot_target_gene(whole.integrated, 'Cd226')
plot_target_gene(whole.integrated, 'Mrc1')

## extract macs

macs <- subset(whole.integrated, idents=c(1:8))
save('macs', file = 'macs/mieloid_custom.RData')

## all genes

plot_markers <- function(object, gene, out_dir, reduction="umap", assay='SCT', slot='data') {
  data <- GetAssayData(object, slot = slot, assay = assay)
  red <- object@reductions[[reduction]]@cell.embeddings
  red <- as.data.frame(red)
  colnames(red) <- paste0(reduction, 1:ncol(red))
  genes_indexes <- which(rownames(data) %in% gene[[1]])
  expression_signaling <- data[genes_indexes,]
  red$expression <- expression_signaling
  onlyBreak <- floor(max(red$expression))
  ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.05)+
    theme_bw(base_size=11) +
    scale_color_gradientn(colours=c("grey", "red", "red3"),breaks=c(0, onlyBreak),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "italic"),
          plot.margin=grid::unit(c(0,0,0.2,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gene)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  ggsave(sprintf("%s/%s.png", out_dir, gene), width = 2.2, height = 2.2, dpi=600, units='in')
}

sapply(top5$gene, function(x) plot_markers(whole.integrated, x, 'review/marker_genes'))
plot_markers(whole.integrated, 'Cxcl14', 'review/marker_genes')
plot_markers(whole.integrated, 'Ptprc', 'review/marker_genes')
plot_markers(whole.integrated, 'Mrc1', 'review/marker_genes')
plot_markers(whole.integrated, 'Treml4', 'review/marker_genes')

sapply(c('S100a8', 'Cd79a', 'Cd3d', 'Nkg7',
         'Lyz2', 'Mrc1', 'Ly6c2', 'Treml4'), function(x) plot_markers(whole.integrated, x, 'review/marker_genes_alex'))


plot_m1 <- function(object, gene, out_dir, reduction="umap", assay='SCT', slot='data') {
  data <- GetAssayData(object, slot = slot, assay = assay)
  red <- object@reductions[[reduction]]@cell.embeddings
  red <- as.data.frame(red)
  colnames(red) <- paste0(reduction, 1:ncol(red))
  genes_indexes <- which(rownames(data) %in% gene[[1]])
  expression_signaling <- data[genes_indexes,]
  red$expression <- expression_signaling
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c('WT', 'KO'))
  red$cluster <- object$integrated_snn_res.0.2
  onlyBreak <- floor(max(red$expression))
  ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.1)+
    theme_bw(base_size=11) +
    facet_grid(.~genotype)+
    scale_color_gradientn(colours=c("grey", "red", "red3"),breaks=c(0, onlyBreak),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    xlim(NA, 5)+
    ylim(0, NA)+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "italic"),
          plot.margin=grid::unit(c(0,0,0.2,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gene)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  ggsave(sprintf("%s/%s.png", out_dir, gene), width = 4.25, height = 2.2, dpi=600, units='in')
}


m1_genes <- c('Cd86', 'Cd80', 'Cd68', 'Tlr2', 'Tlr4', 'Cxcl10')
sapply(m1_genes, function(x) plot_m1(whole.integrated, x, 'review/m1_genes'))


## correct markers


df <- cbind(as.character(whole.integrated$genotype), as.character(whole.integrated$custom_clusters)) %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('genotype', 'cluster'))

df$cluster <- factor(x = df$cluster, levels = as.character(sort(as.numeric(levels(df$cluster)))))
df$genotype <- factor(x = df$genotype, levels = c("WT", "KO"))
p <- df %>% group_by(cluster, genotype) %>% dplyr::summarise(n = n())
per_gt <- table(whole.integrated$genotype) %>% as.data.frame() %>% magrittr::set_colnames(c('genotype', 'count'))

p <- p %>%
  group_by(genotype) %>%
  mutate(total = per_gt$count[match(genotype, per_gt$genotype)]) %>%
  group_by(cluster, add=TRUE) %>%
  mutate(per=round(100*n/total,2))

ggplot(p,aes(x=cluster,y=per, fill=genotype))+ 
  geom_bar(stat = 'identity', position=position_dodge(width = 0.5), width = 0.45) + 
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0, legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  geom_segment(aes(x = 2, y = 8, xend = 2, yend = 6),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 6, y = 13, xend = 6, yend = 11),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 7, y = 13, xend = 7, yend = 11),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 8, y = 10, xend = 8, yend = 8),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 9, y = 10, xend = 9, yend = 8),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  ylab('cells per clutser, %')


ggsave("corrected_hist.png", width = 4, height = 4, dpi=600, units='in')
ggsave("corrected_hist.svg", width = 4, height = 4, dpi=600, units='in')
