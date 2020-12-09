#' Filter genes in RNA-seq ExpressionSet
#'
#' @param eset ExpressionSet with 'counts' assayDataElement and group column in pData
#'
#' @return filtered \code{eset}
#' @export
#'
filter_genes <- function(eset) {
  els <- Biobase::assayDataElementNames(eset)
  els <- ifelse('counts' %in% els, 'counts', 'exprs')
  counts <- Biobase::assayDataElement(eset, els)

  keep <- edgeR::filterByExpr(counts, group = eset$group)
  eset <- eset[keep, ]
  if (!nrow(eset)) stop("No genes with reads after filtering")
  return(eset)
}
