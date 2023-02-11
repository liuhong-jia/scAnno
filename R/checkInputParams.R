##########################################################################################################
#' @title checkPython
#' @return TRUE.

##########################################################################################################

checkPython <- function() {
	py.path <- ifelse(Sys.which("python") == '', readline(prompt = "Please specify python path: "), Sys.which("python"))
	Sys.setenv(RETICULATE_PYTHON = py.path)
	reticulate::use_python(py.path)
	return(TRUE)
}

##########################################################################################################
#' @title checkInputParams
#' @reutur TRUE.

##########################################################################################################

checkInputParams <- function(obj.seu, ref.expr, ref.anno, save.markers, cluster.col) {
	if (!inherits(obj.seu, 'Seurat')) {
		stop(sprintf('%s must be Seurat object, exiting...', deparse(substitute(obj.seu))))
	}
	if (is.null(ref.expr) && is.null(ref.anno)) {
		stop('Cell type and reference expression profile must be provided, exiting...')
	}
	if (!is.null(ref.expr) && is.null(ref.anno)) {
		stop('Cell type must be assigned for refernce profile, exiting...') 
	}
	if (is.null(ref.expr) && !is.null(ref.anno)) {
		stop('Refernce profile must be provided when annotations assigned, exiting...')
	}
	if (!inherits(save.markers, 'character')) {
		stop('Please specified the filename of makers need to be saved, exiting...')
	}
	if (!is.null(ref.expr) && !is.null(ref.anno)) {	
		if (!inherits(ref.expr, 'matrix') && !inherits(ref.expr, 'data.frame')) {
			stop(sprintf('%s must be gene expression matrix or data.frame, exiting...', deparse(substitute(ref.expr))))
		}
		if (!inherits(ref.anno, 'character')) {
			stop(sprintf('%s must be a vector of characters, exiting...', deparse(substitute(ref.anno))))
		}
	}
	if (!(cluster.col %in% colnames(obj.seu@meta.data))) {
		stop(sprintf('%s cannot be detected in the metadata of querySeurat object, exiting...', deparse(substitute(cluster.col))))
	}
	return(TRUE)
}