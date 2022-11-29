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
    
    data(gene.anno)
    
    [INFO] Protein-coding genes.
    
    data(GSE136103)
    
    obj.seu <- GSE136103
    
    [INFO] 35 single-cell clusters that need to be annotated.
    
    ref.obj <- Human_cell_landscape
    
    [INFO] Reference spectrum of human cell landscape, including 63 cell types across more than 50 tissue types.
    
    ref.expr <- GetAssayData(ref.obj, slot = 'data') %>% as.data.frame
    
    ref.anno <- Idents(ref.obj) %>% as.character

# single cell annotation

***

    scAnno <- scAnno (query = obj.seu , ref.expr = ref.expr , ref.anno = ref.anno , gene.anno = gene.anno)
    
    $query
    An object of class Seurat
    21898 features across 16036 samples within 1 assay
    Active assay: RNA (21898 features, 2830 variable features)
    2 dimensional reductions calculated: pca, umap

    $reference
    An object of class Seurat
    17020 features across 5561 samples within 1 assay
    Active assay: RNA (17020 features, 0 variable features)

    $predict_lable
                   C0                    C1                    C2
             "T cell"              "T cell"              "T cell"
                   C3                    C4                    C5
             "T cell"              "T cell"       "M2 Macrophage"
                   C6                    C7                    C8
             "T cell"           "Mast cell"            "CB CD34+"
                   C9                   C10                   C11
         "Neutrophil"          "Basal cell"          "Macrophage"
                  C12                   C13                   C14
             "T cell"    "Endothelial cell"            "Monocyte"
                  C15                   C16                   C17
    "Endothelial cell"    "Endothelial cell"              "T cell"
                  C18                   C19                   C20
         "Macrophage"  "Smooth muscle cell"              "T cell"
                  C21                   C22                   C23
    "Smooth muscle cell"              "B cell"            "Monocyte"
                  C24                   C25                   C26
             "T cell"              "T cell" "B cell (Plasmocyte)"
                  C27                   C28                   C29
        "Dendritic cell"    "Endothelial cell"    "Endothelial cell"
                  C30                   C31                   C32
             "B cell"          "Fibroblast"    "Endothelial cell"
                  C33                   C34
     "Dendritic cell"          "Enterocyte"

    $predict_score
    [1] 0.08963589 0.09507032 0.05481134 0.07018293 0.08915338 0.81135637
    [7] 0.06396457 0.98203585 1.00000612 1.00001372 1.00108067 0.39015408
    [13] 0.07870654 0.17975205 0.22740883 0.11017509 0.25405912 0.09030207
    [19] 0.60888662 0.55502352 0.08643765 0.36434743 0.69280218 0.08491926
    [25] 0.09243409 0.05953276 1.00000000 0.46652109 0.07564509 0.26668928
    [31] 0.30561270 1.00000043 0.11364500 0.25072592 0.99999233

    [INFO] Show annotation results...
![微信图片_20221121111208](https://user-images.githubusercontent.com/115637576/202955413-1362c778-0d2a-419f-9c12-191a63389100.png)

