### Visium data analysis ###

### load packages ###
library(Seurat)
library(ggplot2)

### read in data, create data object, and QC ###
OMRD_3A_assay <- Read10X("./spaceranger_output/OMRD_3A/outs/filtered_feature_bc_matrix")
OMRD_3A <- CreateSeuratObject(counts = OMRD_3A_assay[[1]], project = "OMRD")
OMRD_3A_protein <- CreateAssayObject(counts = OMRD_3A_assay[[2]])
OMRD_3A[['Protein']] <- OMRD_3A_protein
image <- Read10X_Image(image.dir = "./spaceranger_output/OMRD_3A/outs/spatial/")
image <- image[Cells(x = OMRD_3A)]
DefaultAssay(object = image) <- 'RNA'
OMRD_3A[['OMRD_3A']] <- image
OMRD_3A <- RenameCells(OMRD_3A, add.cell.id = "OMRD_3A")
OMRD_3A[["percent.mt"]] <- PercentageFeatureSet(OMRD_3A, pattern = "^MT-")

OMRD_3A_filter <- subset(OMRD_3A, nCount_RNA>=100 & nFeature_RNA>=50 & percent.mt<=25)

### data normalization ###
DefaultAssay(OMRD_3A_filter)<-"RNA"
OMRD_3A_filter <- NormalizeData(OMRD_3A_filter, normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(OMRD_3A_filter)<-"Protein"
OMRD_3A_filter2 <- NormalizeData(OMRD_3A_filter2, normalization.method = "LogNormalize", scale.factor = 100)

### plot features ###
SpatialFeaturePlot(OMRD_3A_filter, features = c("MKI67"), pt.size.factor = 2.7,image.alpha = 0)+ 
  theme(legend.position = "right")


### superpixel level gene expression imputation was performed with iSTAR in python ###

