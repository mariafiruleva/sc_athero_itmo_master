## X_B (add asterics)
library(dplyr)
library(ggplot2)
load("/mnt/tank/scratch/mfiruleva/scn/data/athero_merged/trajectory/without_11/without_4/out/custom_slingshot/rdata/trajectory.RData")
load('/mnt/tank/scratch/mfiruleva/scn/data/athero_merged/trajectory/without_11/out/slingshot/rdata/object.RData')
df <- cbind(as.character(obj$GSE), as.character(obj$old_clusters)) %>% as.data.frame() %>% 
  magrittr::set_colnames(c('study', 'cell_type')) %>% 
  mutate(cell_type = ifelse(cell_type == 0, 'Foam macrophages',
                            ifelse(cell_type == 2, 'Monocytes',
                                   ifelse(cell_type == 5, 'MHC-II macrophages', 'Inflammatory macrophages'))))

df$cell_type <- factor(x = df$cell_type, levels = c('Monocytes', 'Foam macrophages', 'Inflammatory macrophages', 'MHC-II macrophages'))
new_labels <- c(expression('Monocytes', 'Foam macrophages', 'Inflammatory macrophages', ~italic(MHC-II)^{"hi"}~macrophages))
df$study <- factor(x = df$study)


p <- df %>% group_by(cell_type, study) %>% dplyr::summarise(n = n())
per_gt <- table(obj$old_clusters) %>% as.data.frame() %>% magrittr::set_colnames(c('cell_type', 'count')) %>% 
  mutate(cell_type = ifelse(cell_type == 0, 'Foam macrophages',
                            ifelse(cell_type == 2, 'Monocytes',
                                   ifelse(cell_type == 5, 'MHC-II macrophages', 'Inflammatory macrophages'))))
oth_p <- p %>%
  group_by(cell_type) %>%
  mutate(total = per_gt$count[match(cell_type, per_gt$cell_type)]) %>%
  group_by(study, add=TRUE) %>%
  mutate(per=round(100*n/total,2))

ggplot(oth_p,aes(x=cell_type,y=per, fill=study))+ 
  geom_bar(stat = 'identity', position=position_dodge(width = 0.5), width = 0.45) + 
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0,
        legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  ylab('cells: dataset per cluster, %')


ggsave("hist_each_dataset_100perc.png", width = 8, height = 8, dpi=600, units='in')



p <- df %>% group_by(cell_type, study) %>% dplyr::summarise(n = n())
per_gt <- table(obj$GSE) %>% as.data.frame() %>% magrittr::set_colnames(c('study', 'count'))
oth_p <- p %>%
  group_by(study) %>%
  mutate(total = per_gt$count[match(study, per_gt$study)]) %>%
  group_by(cell_type, add=TRUE) %>%
  mutate(per=round(100*n/total,2))
sum(oth_p$per)
sum((oth_p %>% filter(study == 'GSE155513'))$per)
ggplot(oth_p,aes(x=cell_type,y=per, fill=study))+ 
  geom_bar(stat = 'identity', position=position_dodge(width = 0.5), width = 0.45) + 
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0,
        legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  ylab('cells proportions, %')


ggsave("hist_each_dataset_100perc.png", width = 8, height = 8, dpi=600, units='in')

p <- df %>% group_by(cell_type, study) %>% dplyr::summarise(n = n())
per_gt <- table(df$cell_type) %>% as.data.frame() %>% magrittr::set_colnames(c('cell_type', 'count'))

p <- p %>%
  group_by(cell_type) %>%
  mutate(total = per_gt$count[match(cell_type, per_gt$cell_type)]) %>%
  group_by(study, add=TRUE) %>%
  mutate(per=round(100*n/total,2))
sum(p$per)
ggplot(p,aes(x=cell_type,y=per, fill=study))+ 
  geom_bar(stat = 'identity', position=position_dodge(width = 0.5), width = 0.45) + 
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0,
        legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  ylab('cells (dataset per cluster per total no of cells) %')


ggsave("hist_all_clusters_100perc.png", width = 8, height = 8, dpi=600, units='in')

prop_cells <- df %>% group_by(study) %>% dplyr::summarise(n = n()) %>% 
  mutate(percentage = round(100 * (n/ncol(obj)), 2))

ggplot(prop_cells,aes(x=study,y=percentage,fill=study))+ 
  geom_histogram(stat='identity') + 
  geom_text(data=prop_cells, aes(study, percentage+1, label=n), color="black", check_overlap = TRUE)+
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0,
        legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  ylab('cells per myeloid metadataset, %')
ggsave("cells_per_myeloi_meta.png", width = 8, height = 8, dpi=600, units='in')
