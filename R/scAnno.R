##########################################################################################################
#' @title scAnno
#' @param obj.seu Seurat object, which need to be annotated.
#' @param ref.expr Reference gene expression profile. 
#' @param ref.anno Cell type information of reference profile, corresponding to the above `ref.expr`.
#' @param show.plot Show annotated results or not. Default: TRUE.
#' @return Annotated Seurat object.
#' @export 
#' @example	

##########################################################################################################

scAnno <- function(query = obj.seu, ref.expr = ref.expr, ref.anno = ref.anno, show.plot = TRUE, factor.size = 0.1,gene.anno = gene.anno) {
	options(warn = -1)
	python_path <- ifelse(Sys.which("python") == '', readline(prompt = "Please specify python path: "), Sys.which("python"))
	Sys.setenv(RETICULATE_PYTHON = python_path)
	use_python(python_path)
	if (!inherits(obj.seu, 'Seurat')) stop(sprintf('%s must be Seurat object, exiting...', deparse(substitute(obj.seu))))
	if (is.null(ref.expr) && is.null(ref.anno)) stop('Cell type and reference expression profile must be provided, exiting...')
	if (!is.null(ref.expr) && is.null(ref.anno)) stop('Cell type must be assigned for refernce profile, exiting...') 
	if (is.null(ref.expr) && !is.null(ref.anno)) stop('Refernce profile must be provided when annotations assigned, exiting...')
	#gene.anno <- data(gene.anno)
	#readRDS(system.file('external', 'gene.anno.rds', package = 'scAnno'))
	if (!is.null(ref.expr) && !is.null(ref.anno)) {	
	if (!inherits(ref.expr, 'matrix') && !inherits(ref.expr, 'data.frame')) stop(sprintf('%s must be gene expression matrix or data.frame, exiting...', deparse(substitute(ref.expr))))
	if (!inherits(ref.anno, 'character')) stop(sprintf('%s must be a vector of characters, exiting...', deparse(substitute(ref.anno))))
	ref.obj <- CreateSeuratObject(
		counts = ref.expr, 
		project = 'Reference', 
		meta.data = ref.anno %>% as.data.frame(.) %>% `colnames<-`('CellType') %>% `rownames<-`(colnames(ref.expr))
		)
		Idents(ref.obj) <- ref.anno
		ref.obj[["percent.mt"]] <- PercentageFeatureSet(ref.obj, pattern = "^MT-")
		ref.obj[["percent.rb"]] <- PercentageFeatureSet(ref.obj, pattern = "^RP[SL]")
		ref.obj <- subset(
			ref.obj, 
			subset = ( 
				nFeature_RNA > min(quantile(ref.obj@meta.data$nFeature_RNA, 0.05), 300) & 
				nFeature_RNA < max(quantile(ref.obj@meta.data$nFeature_RNA, 0.95), 5000) & 
				nCount_RNA   > min(quantile(ref.obj@meta.data$nCount_RNA, 0.05), 500) & 
				nCount_RNA   < max(quantile(ref.obj@meta.data$nCount_RNA, 0.95), 10000) & 
				percent.mt   < max(quantile(ref.obj@meta.data$percent.mt, 0.95), 15) & 
				percent.rb   < max(quantile(ref.obj@meta.data$percent.rb, 0.95), 35)
			)
		)
	}
	ref.obj <- ref.obj[which(rownames(ref.obj) %in% gene.anno[, 'gene_name']), ]
	avg.expr.ref <- AverageExpression(ref.obj,slot= "data")$RNA
	seed.genes <- identSeedGenes(avg.expr.ref)
	ref.markers <- searchMarkersByCorr(GetAssayData(ref.obj) %>% as.matrix, seed.genes, scale.data = TRUE)
	# Auto annotation using deconvolution and logistic regression methods.
	sig.mat <- buildSigMatrix(ref.markers, avg.expr.ref, min.size = 10, max.size = 200)
	coef.vals <- rlmDewcon(sig.mat, AverageExpression(obj.seu,group.by="seurat_clusters",slot="data")$RNA, weight = NULL)
	if(length(names(table(Idents(ref.obj))))<15){
		signature.gene <- signature_gene(ref.obj,ref.markers)
		lr.pred.res <- lrPredModel(ref.obj[signature.gene, ], obj.seu)
		}
		else{
		markers.lst <- searchMarkersByROC(ref.obj[ref.markers %>% unlist(.) %>% unique,], return.thresh = 0.05)
		lr.pred.res <- lrPredModel(ref.obj[markers.lst %>% unlist(.) %>% unique, ], obj.seu) 
			} 
		lr.pred.score <- lr.pred.res$cellType
		colnames(lr.pred.score) <- paste("X",colnames(lr.pred.score) %>% as.numeric +1,sep="")
		#Merge the scores of the two models
		rownames(coef.vals)<-gsub('_',' ' ,rownames(coef.vals))
		names <-cbind(rownames(coef.vals),names(table(Idents(ref.obj))))
		rownames(lr.pred.score)<-gsub("_|&|\\+|\\.|-|/|\\(|\\)"," ",lr.pred.score %>% as.matrix %>% rownames(.))
		com.names<-intersect(rownames(coef.vals),rownames(lr.pred.score))
		score = -log10(1- coef.vals[com.names,] %>% as.matrix * lr.pred.score[com.names,colnames(coef.vals)]+10^(-10))
		rownames(score) <- names[match(rownames(score),names),2] 
		#rownames(score) <- names(table(Idents(ref.obj)))
		prop <- score/rowSums(score)
		prop[prop<0 | is.na(prop)] <- 0
		rate<-apply(prop,2,max) %>% as.vector  
		lable<-apply(prop,2, function(t) rownames(prop)[which.max(t)]) %>% as.vector
		obj.seu <- AddMetaData(obj.seu,metadata = obj.seu $seurat_clusters %>% as.vector,col.name='scAnno')
		cell.cnts <- table(Idents(obj.seu))
		#Assign new annotation labels to the test set.
		anno.res<-apply(prop,2,which.max) %>% rownames(prop)[.]
		for(i in 0:length(levels(obj.seu$seurat_clusters))-1){
		obj.seu$scAnno[obj.seu$seurat_clusters==i]<-anno.res[i+1]
			}
		if (show.plot) {
		(DimPlot(obj.seu,group.by = 'seurat_clusters',label = TRUE)|
		 DimPlot(obj.seu,group.by = 'scAnno',label = TRUE)) %>% print
			}
		return(list(query = obj.seu, reference = ref.obj, predict_lable = lable , predict_score= rate))
}

