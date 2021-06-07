suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))

set.seed(1)

## PARSE ARGUMENTS

parser <-
  ArgumentParser(description = 'Reanalyze the prepared seurat rdata object')
parser$add_argument('--in_rda',
                    type = "character",
                    help = 'path to the input object')
parser$add_argument('--out_rda',
                    type = "character",
                    help = 'path to the out object')
parser$add_argument('--clusters',
                    type = "integer", nargs='+',
                    help = 'target clusters')
parser$add_argument('--ident',
                    type = "character",
                    help = 'resolution name')
parser$add_argument('--annot',
                    type = "character",
                    help = 'path to the annotation file')

args <- parser$parse_args()
ident <- args$ident

print(args)

## LOAD OBJECT, SET IDENT, SUBSET

load(args$in_rda)

Idents(whole.integrated) <- args$ident

obj <- subset(whole.integrated, idents = args$clusters)

if (grepl('GSE131776', args$in_rda)) {
  sub_cells <- Embeddings(whole.integrated, reduction = 'umap') %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'cells') %>% 
    filter(UMAP_1 > 6.978 & UMAP_1 < 8.14 & UMAP_2 > 3.036 & UMAP_2 < 4.606)
  obj <- subset(obj, cells=sub_cells$cells, invert=T)
}

## ADD META DATA

annot <- fread(args$annot)

for (col in colnames(annot)[!grepl('SRS', colnames(annot))]) {
  obj[[col]] <- annot[[col]][match(obj$sample, annot$SRS)]
}
obj[['organism']] <- 'Mus musculus'

## SAVE

save(obj, file=args$out_rda)
