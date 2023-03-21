# scAnno(**sc**RNA-seq data **A**nnotation)

***

scAnno is an automated annotation tool for single-cell RNA sequencing datasets primarily based on the single cell cluster levels, using a joint deconvolution strategy and logistic regression. We explicitly created a complete reference atlas (reference expression profiles) of 30 cell types from the Human Cell Landscape (HCL) covering more than 50 human tissues to support this novel methodology (scAnno). scAnno offers a possibility to obtain genes with high expression and specificity in a given cell type as cell type-specific genes (marker genes) by combining co-expression genes with seed genes as a core. Of importance, scAnno can accurately identify cell type-specific genes based on cell type reference expression profiles without prior information. 
![1679397927741](https://user-images.githubusercontent.com/115637576/226592587-128e2d85-ffc6-4fce-ad09-830a14c424f9.png)

# Installing the package

***
 
To install scAnno,we recommed using devtools:  

    #install.packages("devtools")  
    devtools::install_github("liuhong-jia/scAnno")  

***

# Dependencies
- R version >= 3.5.0.
- R packages: Seurat, dplyr, reticulate, MASS, irlba, future, progress, parallel, glmnet, knitr, rmarkdown, devtools

# Guided Tutorials
For this tutorial, we apply the human single cell reference atlas(Human_cell_landscape) built into the package to predict a scRNA-seq dataset(GSE136103) derived from liver tissue that has been processed by the standard Seurat process and entered as a query object.


# Prepare input data

    library(scAnno)
    
    #Import human cell type reference profile.
    data(Human_cell_landscape)
    
    #Import protein coding gene(19814 genes) to filter reference expression profile.
    data(gene.anno)
    
    #Import TCGA bulk data in pan-cancer.
    data(tcga.data.u)
    
    #A liver tissue data set to be annotated.
    data(GSE136103)
    
  
# Set parameters
|**Parameters**|**Description**                      |
|----------|-----------------------------------------|
|query     |Seurat object, which need to be annotated|
|ref.expr  |Reference gene expression profile.       |
|ref.anno  |Cell type information of reference profile, corresponding to the above `ref.expr`.|
|save.markers|Specified the filename of makers need to be saved.Default: markers.|
|cluster.col|Column name of clusters to be annotated in meta.data slot of query Seurat object. Default: seurat_clusters.|
|factor.size|Factor size for scaling the weight of gene expression. Default: 0.1.|
|seed.num|Number of seed genes of each cell type for recognizing candidate markers, only used when method = 'co.exp'. Default: 10.|
|redo.markers|Re-search candidate markers or not. Default: FALSE.|
|gene.anno|Gene annotation data.frame. Default: gene.anno.|
|permut.num|Number of permutations for estimating p-values of annotations. Default: 100.|
|show.plot|Show annotated results or not. Default: TRUE.|
|verbose|Show running messages or not. Default: TRUE.|
|tcga.data.u|bulk RNA-seq data of pan-cancer in TCGA.|

**Note**: The parameter save.markers means that the marker genes will be stored in a temporary file, so that the next time the same reference expression is used, it will not have to be run again.
# single cell annotation

    # Seurat object, which need to be annotated.
    obj.seu <- GSE136103
    
    #Seurat object of reference gene expression profile.
    ref.obj <- Human_cell_landscape
    
    #Reference gene expression profile.
    ref.expr <- GetAssayData(ref.obj, slot = 'data') %>% as.data.frame
    
    #Cell type information of reference profile, corresponding to the above `ref.expr`.
    ref.anno <- Idents(ref.obj) %>% as.character
    
# Run scAnno to annotate cell 
Details of the results is described in the table below.
|**output**|**details**|
|------|-------|
|query|Seurat object, which need to be annotated.|
|reference|Seurat object of reference gene expression profile.|
|pred.label|Cell types corresponding to each cluster.|
|pred.score|The prediction score for each cluster,corresponding to `pred.label`.|
|pvals|Significance level of the predicted scores, corresponding to `pred.score`.|

	results = scAnno(query = obj.seu,
		ref.expr = ref.expr,
		ref.anno = ref.anno,
		save.markers = "markers",
		cluster.col = "seurat_clusters",
		factor.size = 0.1,
		pvalue.cut = 0.01,
		seed.num = 10,
		redo.markers = FALSE,
		gene.anno = gene.anno,
		permut.num = 100,
		show.plot = TRUE,
		verbose = TRUE,
		tcga.data.u = tcga.data.u
		)
	[INFO] Checking the legality of parameters
	[INFO] 30 cell types in reference, 35 clusters in query objects
	[INFO] Deconvolution by using RLM method
	[INFO] Logistic regression for cell-type predictions, waiting...
	[INFO] Merging the scores of both models, and assign annotations to clusters
	[INFO] Estimating p-values for annotations...
	[INFO] Finish!
	
	results$query
	An object of class Seurat
	21898 features across 16036 samples within 1 assay
	Active assay: RNA (21898 features, 2830 variable features)
	2 dimensional reductions calculated: pca, umap
	
	results$reference
	An object of class Seurat
	17020 features across 5561 samples within 1 assay
	Active assay: RNA (17020 features, 0 variable features)

    results$pred.label
              C0                    C1                    C2
             "T cell"              "T cell"              "T cell"
                   C3                    C4                    C5
             "T cell"              "T cell"      "Dendritic cell"
                   C6                    C7                    C8
             "T cell"              "T cell"              "T cell"
                   C9                   C10                   C11
           "Monocyte"     "Epithelial cell"          "Macrophage"
                  C12                   C13                   C14
             "T cell"    "Endothelial cell"            "Monocyte"
                  C15                   C16                   C17
   		"Endothelial cell"  "Endothelial cell"       "T cell"
                  C18                   C19                   C20
         "Macrophage"    "Smooth muscle cell"            "T cell"
                  C21                   C22                   C23
 	    "Smooth muscle cell"       "B cell"            "Monocyte"
                  C24                   C25                   C26
             "T cell"             "T cell"  "B cell (Plasmocyte)"
                  C27                   C28                   C29
     	"Dendritic cell"    "Endothelial cell"  "Endothelial cell"
                  C30                   C31                   C32
              "B cell"        "Stromal cell"    "Endothelial cell"
                  C33                   C34
      "Dendritic cell"     "Epithelial cell"
	
    results$pred.score
	[1] 0.9999944 0.9987957 0.9972636 1.0000000 1.0000000 0.9979784 0.9999558
	[8] 0.9977544 0.9981475 0.9995275 0.8849131 0.9949165 0.9999916 0.9986469
	[15] 0.9586247 0.9774114 0.9978670 1.0000000 0.9902297 0.9995580 1.0000000
	[22] 0.9999576 0.9979638 0.6009568 1.0000000 0.9984880 0.9985971 1.0000000
	[29] 0.9984917 0.9989927 0.9806239 0.6151706 0.9796672 1.0000000 1.0000000

    
    results$pvals
     	C0           C1           C2           C3           C4           C5
	4.007509e-17 2.970442e-16 3.538215e-15 3.969419e-17 3.969419e-17 0.000000e+00
          C6           C7           C8           C9          C10          C11
	4.278350e-17 1.616284e-15 8.570515e-16 0.000000e+00 0.000000e+00 0.000000e+00
         C12          C13          C14          C15          C16          C17
	4.026227e-17 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 3.969419e-17
         C18          C19          C20          C21          C22          C23
	0.000000e+00 0.000000e+00 3.969419e-17 0.000000e+00 0.000000e+00 0.000000e+00
         C24          C25          C26          C27          C28          C29
	3.969419e-17 4.921965e-16 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
         C30          C31          C32          C33          C34
	0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00

# Visualization
Show annotation results...The left graph represents the UMAP plot of cluster of query dataset，and the right graph represents the annotation of scAnno.

    DimPlot(results$query, group.by = "seurat_clusters", label = TRUE) | DimPlot(results$query, group.by = 'scAnno', label = TRUE)
    
![c75c451c3a993d5f5a78adae32947c4](https://user-images.githubusercontent.com/115637576/218242912-44df6b81-7501-4840-aa1d-d97bb7121aea.png)

