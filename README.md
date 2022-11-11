# scAnno

***

A deconvolution strategy-based automatic cell type annotation tool for single-cell RNA sequencing datasets.
![image](https://user-images.githubusercontent.com/115637576/199387392-3d66fb26-e3d3-43c9-9378-04f541600e3f.png)


# Installing the package

***

You can install the package using devtools::install_github:

    devtools::install_github("liuhong-jia/scAnno")

# Getting started with scAnno

***

In this tutorial we will use GSE136103 (Liver) as an example.

    library(scAnno)
    library(Seurat)
    library(patchwork)
    library(dplyr)
    library(MASS)
    library(irlba)
    library(future)
    library(glmnet)
    library(reticulate)
    library(parallel)
    
    data(Human_cell_landscape)
    
    data(GSE136103)
    
    obj.seu <- GSE136103
    
    [INFO] 35 single-cell clusters that need to be annotated.
    
    ref.obj <-Human_cell_landscape
    
    [INFO] Reference spectrum of human cell landscape, including 63 cell types across more than 50 tissue types.
    
    ref.expr <- GetAssayData(ref.obj, slot = 'data') %>% as.data.frame
    
    ref.anno <- Idents(ref.obj) %>% as.character

# single cell annotation

***

    scAnno <- scAnno (query = obj.seu , ref.expr = ref.expr , ref.anno = ref.anno , gene.anno = gene.anno)
    
    scAnno
    $query
    An object of class Seurat
    21898 features across 16036 samples within 1 assay
    Active assay: RNA (21898 features, 2830 variable features)
    2 dimensional reductions calculated: pca, umap
    
    $reference
    An object of class Seurat
    17020 features across 11981 samples within 1 assay
    Active assay: RNA (17020 features, 0 variable features)
    
    $predict_lable
    [1] "T cell"                              "T cell"
    [3] "T cell"                              "T cell"
    [5] "T cell"                              "M2 Macrophage"
    [7] "T cell"                              "Mast cell"
    [9] "Erythroid progenitor cell (RP high)" "Neutrophil"
    [11] "Enterocyte"                          "Macrophage"
    [13] "T cell"                              "Endothelial cell (APC)"
    [15] "Antigen presenting cell (RPS high)"  "Sinusoidal endothelial cell"
    [17] "Endothelial cell"                    "T cell"
    [19] "Macrophage"                          "Smooth muscle cell"
    [21] "T cell"                              "Smooth muscle cell"
    [23] "B cell"                              "Monocyte"
    [25] "T cell"                              "Mast cell"
    [27] "B cell (Plasmocyte)"                 "Dendritic cell"
    [29] "Endothelial cell (APC)"              "Sinusoidal endothelial cell"
    [31] "B cell"                              "Mesothelial cell"
    [33] "Endothelial cell"                    "Dendritic cell"
    [35] "Sinusoidal endothelial cell"
    
    $predict_score
    [1] 0.07787374 0.07216236 0.04723584 0.08651074 0.09659525 0.70952887
    [7] 0.08282688 0.86187648 1.00015873 1.00004725 1.00128895 0.45628805
    [13] 0.08370610 0.25861030 0.66387577 0.35116881 0.33880931 0.07431788
    [19] 0.54281711 0.58160045 0.09688529 0.32773497 0.67251650 0.34593674
    [25] 0.07810757 0.13176582 1.00000001 0.36813368 0.27101113 0.33095474
    [31] 0.32572143 1.00001112 0.38655282 0.30972457 0.17550569

   


    [INFO] Show annotation results...
![2e84df1c21a02be5a72e2768c5580ae](https://user-images.githubusercontent.com/115637576/201280348-864beccb-3797-46b4-aa20-4305ec1a9416.png)
