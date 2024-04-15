
# This script is for better visualiztion of dot plots after clustering of Seurat objects' data

# It inspired from Ming Tang and the Eye Bioinformatician

# Load libraries

library(Seurat)
library(tidyverse)
library(presto)
library(ComplexHeatmap)
library(circlize)



# So.. to start the visualization we do the basic initial analysis (standard processing) on scRNAseq data -->

# 1. Create Seurat object
# 2. Deal with mt genes
# 3. Normalize data
# 4. Find the variable features
# 4. Scale data
# 5. Run PCA 
# 6. Find neighbours
# 7. Find clusters
# 8. Run UMAP and START TO VISUALIZE YOUR DATA
# 9. Also find marker genes for each cluster

# For example

markers = presto::wilcoxauc(myObject, "seurat_clusters", assay = "data")
markers = top_markers(markers, n = 10, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
markers
# These marker genes can be visualized in dot plots
all_markers =  markers %>%
               select(-rank) %>%
               unclass() %>%
               stack() %>%
               pull(values) %>%
               unique() %>%
               .[!is.na(.)]
# In Seurat the visualization is like this:

p = Dotplot( object = myObgect, feauture = all_markers)


# To start reproducing the analysis and visualization

df = p$data
head(df)

# Create a matrix for the scaled expression 
exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

head(exp_mat)

# Create a matrix of the cells expressing a gene

percent_mat = df %>%
              select(-avg.exp, -avg.exp.scaled) %>%
              pivot_wider(names_from = id, values_from = prct.exp) %>%
              as.data.frame()
 row.names(percent_mat) = percent_mat$features.plot
 percent_mat =  percent_mat[,-1] %>% as.matrix()
 
 head(percent_mat)
 
# Check that the range is from 0 - 100
 range(percent_mat)
 
# Check that the two matrices have the same nnumber of dimensions
 dim(exp_mat)
 dim(percent_mat)
 
 
 # Load libraries t omakes things more colourful
 library(viridis)
 #INSTALL Polychrome library
 install.packages("Polychrome", repos="http://R-Forge.R-project.org")
 library(Polychrome)
 Polychrome::swatch(viridis(40))
 
 ## get an idea of the ranges of the matrix
 quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))
 
 ## any value that is greater than 2 will be mapped to yellow
 col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])
 
 
 cell_fun = function(j, i, x, y, w, h, fill){
   grid.rect(x = x, y = y, width = w, height = h, 
             gp = gpar(col = NA, fill = NA))
   grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
               gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
 
 ## also do a kmeans clustering for the genes with k = 4
 Heatmap(exp_mat,
         heatmap_legend_param=list(title="expression"),
         column_title = "clustered dotplot", 
         col=col_fun,
         rect_gp = gpar(type = "none"),
         cell_fun = cell_fun,
         row_names_gp = gpar(fontsize = 5),
         row_km = 4,
         border = "black")
 
# One can now take advantage of all the awesome annotation function 
# to add annotation bars in complexheatmap. 
 colnames(exp_mat)
 
 
library(RColorBrewer)
cluster_anno<-  c("CD4T", "B", "CD4T", "Mono", "NK", "CD8T", "CD14_Mono", "DC", "Platelet")
 
 column_ha<- HeatmapAnnotation(
   cluster_anno = cluster_anno,
   col = list(cluster_anno = setNames(brewer.pal(8, "Paired"), unique(cluster_anno))
   ),
   na_col = "grey"
 )
 
Heatmap(exp_mat,
         heatmap_legend_param=list(title="expression"),
         column_title = "clustered dotplot", 
         col=col_fun,
         rect_gp = gpar(type = "none"),
         cell_fun = cell_fun,
         row_names_gp = gpar(fontsize = 5),
         row_km = 4,
         border = "black",
         top_annotation = column_ha)

# use grid.circle() in both the heatmap body and the legend

## To make the size of the dot in the heatmap body comparable to the legend, I used fixed
## size unit.(2, "mm) rather than min(unit.c(w, h).

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= pindex(percent_mat, i, j)/100 * unit(2, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}


lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  ))

set.seed(123)    
hp<- Heatmap(exp_mat,
             heatmap_legend_param=list(title="expression"),
             column_title = "clustered dotplot", 
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 5),
             row_km = 4,
             border = "black",
             top_annotation = column_ha)

draw( hp, annotation_legend_list = lgd_list)


# Note that, we use the percentage as the radius of the circle. 
# The area of the dot is pi *r^2. 
# It will appear 4 times different if the original percentage is 0.5 versus 1 (2 times difference) which may not be desirable.

# We can do a square root of the radius:

layer_fun1 = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= sqrt(pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}

lgd_list1 = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0  * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  ))

set.seed(123)    
hp1<- Heatmap(exp_mat,
              heatmap_legend_param=list(title="expression"),
              column_title = "clustered dotplot", 
              col=col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun1,
              row_names_gp = gpar(fontsize = 5),
              row_km = 4,
              border = "black",
              top_annotation = column_ha)

draw( hp1, annotation_legend_list = lgd_list1)

# use grid.points() in both the heatmap body and the legend

## note for grid.points, use col for the filling color of the points, while in grid.circle, use fill for the filling color of the circle. I should learn more about {grid}

layer_fun2 = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(exp_mat, i, j))),
              size = pindex(percent_mat, i, j)/100 * unit(4, "mm"),
              pch = 16
  )
}


lgd_list2 = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

set.seed(123)
hp2<- Heatmap(exp_mat,
              heatmap_legend_param=list(title="expression"),
              column_title = "clustered dotplot", 
              col=col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun2,
              row_names_gp = gpar(fontsize = 5),
              row_km = 4,
              border = "black",
              top_annotation = column_ha)

draw( hp2, annotation_legend_list = lgd_list2)



