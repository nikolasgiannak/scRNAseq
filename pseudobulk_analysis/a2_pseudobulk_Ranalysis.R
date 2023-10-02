
################################################################################
######################### Load libraries #######################################
################################################################################
library(here)
library(stringr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(purrr)
library(readr)
library(harmony)
library(scCustomize)
library(SeuratDisk)
library(ggpointdensity)
library(viridis)
library(grid)

################################################################################
############ read in counts files fromm the a1_Get_the_pub_avail_data.sh #######
################################################################################

## the count matrix is cell x gene
read_counts<- function(file){
  x<- read_csv(file)
  x<- as.data.frame(x)
  cells<- x$index
  
  mat<- as.matrix(x[,-1])
  
  rownames(mat)<- cells
  mat<- t(mat)
  return(mat)
}


counts_files<- list.files("~/blog_data/GSE154763", full.names = TRUE, pattern = "*raw_count.txt")
counts_files

# only use the last 3 cancer types to demonstrate, as the data are too big. 
counts_files<- counts_files %>%
  tail(n=3)

samples<- map_chr(counts_files, basename) 
samples<- str_replace(samples, "(GSE[0-9]+)_(.+)_raw_count.txt", "\\2")
names(counts_files)<- samples
counts<- purrr::map(counts_files, read_counts)

################################################################################
#########read in metadata for each cancer#######################################
################################################################################

## the tissue column contain some non-tumor cells
read_meta<- function(file){
  y<- read_csv(file, guess_max = 100000, 
               col_types = cols(tissue = col_character()))
  y<- as.data.frame(y)
  cells<- y$index
  y<- y[,-1]
  rownames(y)<- cells
  return(y)
}


meta_files<- list.files("~/blog_data/GSE154763", full.names = TRUE, pattern = "metadata.csv")
# take only 3 cancer types
meta_files<- meta_files %>% tail(n=3)
meta_names<- map_chr(meta_files, basename)
meta_names<- str_replace(meta_names, "(GSE[0-9]+)_(.+)_metadata.csv", "\\2")
names(meta_files)<- meta_names

meta<- purrr::map(meta_files, read_meta)

################################################################################
################## Different cancer types have different number of genes. ######
################################################################################

map(counts, nrow)

################################################################################
#################### How many genes are in common among the 3 datasets? ########
################################################################################

map(counts, ~rownames(.x)) %>%
  purrr::reduce( function(x,y) {intersect(x,y)}) %>%
  length()

################################################################################
#####################  make Seurat objects #####################################
################################################################################

objs<- purrr::map2(counts, meta,  
                   ~ CreateSeuratObject(counts = as(.x, "sparseMatrix"), 
                                        meta.data = .y))

## remain only the tumor samples with tissue == T
#subset(objs$ESCA, subset = tissue == "T")

objs<- purrr::map(objs, ~ subset(.x, subset = tissue == "T"))

objs

################################################################################
#############  It is a list of 3 seurat objects. ###############################
################################################################################

preprocessSeurat<- function(obj){
  obj<-
    NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures( selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = "patient", dims.use = 1:30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.6)
  
  return(obj)
}


objs<- purrr::map(objs, preprocessSeurat)


# visualization of Dimplot 

p1<- scCustomize::DimPlot_scCustom(objs$UCEC)

p2<- scCustomize::DimPlot_scCustom(objs$UCEC, group.by = "MajorCluster")

p1


################################################################################

################# Create psedobulk for the single cell data ####################

################################################################################




################################################################################
#### Clean up the cell type annotation so it is consistent across cancer type.###
################################################################################

purrr::map(objs, ~.x$MajorCluster %>% unique() %>% sort())

clean_cluster_name<- function(obj){
  annotation<- obj@meta.data %>%
    mutate(annotation = str_replace(MajorCluster, "^M[0-9]{2}_", "")) %>%
    pull(annotation)
  obj$annotation<- annotation
  return(obj)
}


objs<- purrr::map(objs, clean_cluster_name)

purrr::map(objs, ~.x$annotation %>% unique() %>% sort())
 

################################################################################
#       use presto for more low-level analysis instead of
#    using the wrapper functions in muscat. To use muscat,
# you will need to convert the Seurat object to a SingleCellExperiment object.
################################################################################
library(presto)

objs$PAAD@meta.data %>% 
  head()

################################################################################
#   We will need to create a pseduobulk per cancer type, per sample,
#     and per cell type. those are in the metadata column:
#       ‘annotation’, ‘patient’, ‘cancer’.
################################################################################

create_pseduo_bulk<- function(obj){
  data_collapsed <- presto::collapse_counts(obj@assays$RNA@counts, 
                                            obj@meta.data, 
                                            c('annotation', 'patient', 'cancer'))
  meta_data<- data_collapsed$meta_data 
  
  mat<- data_collapsed$counts_mat
  
  colnames(mat)<- paste(meta_data$cancer, colnames(mat), sep="_")
  return(list(mat = mat, meta_data = meta_data))
}

pseduo_bulk_objs<- purrr::map(objs, create_pseduo_bulk)

### different number of genes per dataset 
purrr::map(pseduo_bulk_objs, ~nrow(.x$mat)) 


genes_per_data<- purrr::map(pseduo_bulk_objs, ~rownames(.x$mat)) 

## common number of genes
common_genes<- purrr::reduce(genes_per_data, intersect)

length(common_genes)


subset_common_mat<- function(mat){
  mat<- mat[common_genes, ]
  return(mat)
}


mats<- map(pseduo_bulk_objs, ~subset_common_mat(.x$mat))
map(mats, dim)

meta_data<- map(pseduo_bulk_objs, ~.x$meta_data)

meta_data$PAAD %>% head()

mats$PAAD[1:5, 1:6]


################################################################################
########################   PCA analysis  #######################################
################################################################################

#   merge all the count table into a big one
mats<- purrr::reduce(mats, cbind)
meta_data<- purrr::reduce(meta_data, rbind)
dim(mats)


final_meta<- meta_data

## log normalize the count matrix 
total_reads<- colSums(mats)

final_mat<- t(t(mats)/total_reads)
final_mat<- log2(final_mat + 1)

library(genefilter)

# choose the top 1000 most variabel genes 
top_genes<- genefilter::rowVars(final_mat) %>% 
  sort(decreasing = TRUE) %>%
  names() %>%
  head(1000)

# subset only the top 1000 genes
expression_mat_sub<- final_mat[top_genes, ]


# calculate the PCA
pca<- prcomp(t(expression_mat_sub),center = TRUE, scale. = TRUE) 

# check the order of the samples are the same.
all.equal(rownames(pca$x), colnames(expression_mat_sub))

PC1_and_PC2<- data.frame(PC1=pca$x[,1], PC2= pca$x[,2])

PC1_and_PC2<- cbind(PC1_and_PC2, final_meta)

head(PC1_and_PC2)

# plot PCA
library(ggplot2)

p1<- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = cancer)) +
  theme_bw(base_size = 14) 

# color it by cell type

library(Polychrome)
set.seed(123)

length(unique(PC1_and_PC2$annotation))


mypal <- kelly.colors(14)[-1]

p2<- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = annotation)) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = mypal %>% unname()) 

p1
p2

# Visualize in heatmap.
library(RColorBrewer)
library(ComplexHeatmap)
#RColorBrewer::display.brewer.all()
set.seed(123)
cols<- RColorBrewer::brewer.pal(n = 3, name = "Set1")
unique(final_meta$cancer)

rownames(final_meta)<- colnames(expression_mat_sub)

ha<- HeatmapAnnotation(df = final_meta[, -2],
                       col=list(cancer = setNames(cols, unique(final_meta$cancer) %>% sort()),
                                annotation= setNames(unname(mypal), unique(final_meta$annotation) %>% sort())),
                       show_annotation_name = TRUE)

col_fun<- circlize::colorRamp2(c(-2,0, 2), colors = c("blue", "white", "red"))

Heatmap(t(scale(t(expression_mat_sub))), top_annotation = ha,
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        name = "expression",
        column_dend_reorder = TRUE,
        raster_quality = 3,
        use_raster = FALSE,
)


# Let’s run Harmony to correct the batches of different cancer type and patients.
set.seed(123)
harmony_embeddings <- harmony::HarmonyMatrix(
  expression_mat_sub, final_meta, c('cancer', 'patient'), do_pca = TRUE, verbose= TRUE
)

rownames(harmony_embeddings)<- rownames(final_meta)
harmony_pc<- data.frame(harmony1=harmony_embeddings[,1], harmony2= harmony_embeddings[,2])

harmony_pc<- cbind(harmony_pc, final_meta)

## plot PCA plot
library(ggplot2)

p1<- ggplot(harmony_pc, aes(x=harmony1, y=harmony2)) + 
  geom_point(aes(color = cancer)) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = cols)

p2<- ggplot(harmony_pc, aes(x=harmony1, y=harmony2)) + 
  geom_point(aes(color = annotation)) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = mypal %>% unname())

p1
p2

# Let’s visulize the heatmap by the sample correlations in the harmony space.

Heatmap(cor(t(harmony_embeddings)), top_annotation = ha, show_row_names = FALSE,
        show_column_names = FALSE, column_dend_reorder = TRUE, 
        name="harmony space\ncorrelation",
        raster_quality = 3,
        use_raster = TRUE)

install.packages("harmony")
devtools::install_github("immunogenomics/harmony", build_vignettes=TRUE)



