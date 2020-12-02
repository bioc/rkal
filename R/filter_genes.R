#' Filter genes in RNA-seq ExpressionSet
#'
#' @param eset ExpressionSet with 'counts' assayDataElement and group column in pData
#'
#' @return filtered \code{eset}
#' @export
#' @keywords internal
#'
filter_genes <- function(eset) {
  counts <- Biobase::assayDataElement(eset, 'counts')
  if (is.null(counts)) return(eset)
  keep <- edgeR::filterByExpr(counts, group = eset$group)
  eset <- eset[keep, ]
  if (!nrow(eset)) stop("No genes with reads after filtering")
  return(eset)
}
