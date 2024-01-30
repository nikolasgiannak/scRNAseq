

if (!requireNamespace("remotes")) install.packages("remotes")
library(remotes)


# Install STACAS for semi-supervised integration of scRNAseq data
remotes::install_github("carmonalab/STACAS")
library(STACAS)

# Load needed libraries and set.seed
install.packages("SeuratObject")
library(SeuratObject)
library(Seurat)
library(dplyr)
library(ggplot2)
library(STACAS)

seed = 1234
set.seed(seed)

if (!requireNamespace("SeuratData")) install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
install.packages("scCustomize")

library(scCustomize)

InstallData("pbmcsca")
data("pbmcsca")

pmbcsca.updated = UpdateSeuratObject(object =pbmcsca)


# Integrate scRNA-seq datasets generated in different batches (in this example, using different methods/technologies)
pbmcsca.integrated <- NormalizeData(pmbcsca.updated) |>
  SplitObject(split.by = "Method")|>
  Run.STACAS()

pbmcsca.integrated <- RunUMAP(pbmcsca.integrated, dms = 1:30) 

# Visualize
Dimplot_1_method = DimPlot(pbmcsca.integrated, group.by = c("Method")) 

Dimplot_1_CellType = DimPlot(pbmcsca.integrated, group.by = c("CellType"), label.size = 0.5) 

DimPlot(pbmcsca.integrated, group.by = c("Experiment")) 

Dimplot_1_Cluster = DimPlot(pbmcsca.integrated, group.by = c("Cluster")) 





# Semi-supervised integration
pbmcsca.semisup <- NormalizeData(pmbcsca.updated) |>
  SplitObject(split.by = "Method")|>
  Run.STACAS(cell.labels = "CellType")

pbmcsca.semisup <- RunUMAP(pbmcsca.semisup, dims = 1:30) 


# Visualize
Dimplot_2_method = DimPlot(pbmcsca.semisup, group.by = c("Method"), label = FALSE) & NoLegend()

Dimplot_2_CellType = DimPlot_scCustom(pbmcsca.semisup, group.by = c("CellType"),
                                      pt.size = 0.005,
                                      colors_use =  DiscretePalette_scCustomize(num_colors = 15, palette = "ditto_seq" )
                                      ) & NoLegend()

DimPlot(pbmcsca.semisup, group.by = c("Experiment")) 

Dimplot_2_Cluster = DimPlot_scCustom(pbmcsca.semisup, 
                            group.by = c("Cluster"),
                            label.size = 0.1,
                            colors_use =  DiscretePalette_scCustomize(num_colors = 15, palette = "varibow" )
                            ) & NoLegend()

Dimplot_2_Cluster = DimPlot_scCustom(seurat_object = pbmcsca.semisup, 
                                     colors_use =  DiscretePalette_scCustomize(num_colors = 15, palette = "ditto_seq" ))


 library(patchwork)
Dimplot_1_method | Dimplot_2_method
Dimplot_1_CellType | Dimplot_2_CellType
Dimplot_1_Cluster | Dimplot_2_Cluster

 (Dimplot_1_method | Dimplot_2_method) / (Dimplot_1_CellType | Dimplot_2_CellType)
