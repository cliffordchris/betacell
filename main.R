#If needed to download libraries "SeuratDisk" and "SeuratData"
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

devtools::install_github('satijalab/seurat-data')

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SeuratData)

  # For local files:
  
  expression_matrix <- ReadMtx(
    mtx = "endocrine.counts.mtx.gz", features = "features.tsv.gz",
    cells = "barcodes.tsv.gz"
  )
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  #can uncomment below code to import saved data
  #seurat_object <- LoadH5Seurat("SeuratProject.h5Seurat")
  
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") 
  #save the current seurat object in the working directory to prevent future imports
  SaveH5Seurat(seurat_object)
  
  #print general characteristics
  seurat_object
  
  # Visualize QC metrics as a violin plot, stores as a png file named "QC plot.png"
  png(file = "QC plot.png", width = 2000, height = 800)
  VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  # # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  # 
  # plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  # plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # CombinePlots(plots = list(plot1, plot2))
  # 
  # #scaling the data
  # all.genes <- rownames(seurat_object)
  # seurat_object <- ScaleData(seurat_object, features = all.genes)
  # #save the current seurat object in the working directory to prevent future imports
  # SaveH5Seurat(seurat_object, filename = "post-scaling") 
  # #Perform linear dimensional reduction
  # seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  # ElbowPlot(seurat_object)
  # 