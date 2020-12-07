.datatable.aware = TRUE

#' Load bulk RNA-Seq data into an ExpressionSet.
#'
#' @param data_dir Directory with raw and quantified RNA-Seq files.
#' @param species Character vector indicating species. Genus and species should be space seperated, not underscore. Default is \code{Homo sapiens}.
#' @param release EnsemblDB release. Should be same as used in \code{\link[drugseqr.data]{build_kallisto_index}}.
#' @param load_saved If TRUE (default) and a saved \code{ExpressionSet} exists, will load from disk.
#' @param save_eset If TRUE (default) and either \code{load_saved} is \code{FALSE} or a saved \code{ExpressionSet} does not exist,
#'   then an ExpressionSet will be saved to disk. Will overwrite if already exists.
#'
#' @return \code{\link[Biobase]{ExpressionSet}} with attributes/accessors:
#' \itemize{
#'   \item \code{sampleNames} from names of raw RNA-Seq files (excluding .fastq.gz suffix).
#'   \item \code{annotation} Character vector of annotation package used.
#'   \item \code{exprs} Length scaled counts generated from abundances for use in
#'     \code{\link[limma]{voom}} (see \code{vignette("tximport", package = "tximport")}).
#'   \item \code{abundance, counts, length} accessed e.g. by \code{assayDataElement(eset, 'length')}.
#'   Imported for exploratory data analysis with DESeq2 variance stabilization transforms.
#'   \item \code{phenoData} added columns:
#'     \itemize{
#'       \item \code{lib.size} library size from \code{\link[edgeR]{calcNormFactors}}.
#'       \item \code{norm.factors} library normalization factors from \code{\link[edgeR]{calcNormFactors}}.
#'     }
#' }
#'
#' @export
#' @seealso \link{get_vsd}
#'
load_seq <- function(data_dir, species = 'Homo sapiens', release = '94', load_saved = TRUE, save_eset = TRUE) {

  # check if already have
  eset_path  <- file.path(data_dir, 'eset.rds')
  if (load_saved & file.exists(eset_path))
    return(readRDS(eset_path))

  # import quants and filter low counts
  q <- import_quants(data_dir, species = species, release = release)

  # construct eset
  annot <- get_ensdb_package(species, release)
  fdata <- setup_fdata(species, release)
  eset <- construct_eset(q$quants, fdata, annot, q$txi.deseq)

  # save eset and return
  if (save_eset) saveRDS(eset, eset_path)
  return(eset)
}

#' Get variance stabilized data for exploratory data analysis
#'
#' @param eset ExpressionSet loaded with \link{load_seq}.
#'   Requires group column in \code{pData(eset)} specifying sample groupings.
#' @param rlog_cutoff Sample number above which will use \code{\link[DESeq2]{vst}}.
#'   instead of \code{\link[DESeq2]{rlog}}. Default is 50.
#'
#' @return \code{DESeqTransform} with variance stabilized expression data.
#' @export
get_vsd <- function(eset, rlog_cutoff = 50) {

  trans_fun <- ifelse(ncol(eset) > rlog_cutoff, DESeq2::vst, DESeq2::rlog)
  pdata <- Biobase::pData(eset)
  txi.deseq <- list(countsFromAbundance = 'no',
                    abundance = Biobase::assayDataElement(eset, 'abundance'),
                    counts = Biobase::assayDataElement(eset, 'counts'),
                    length = Biobase::assayDataElement(eset, 'length')
  )

  dds <- DESeq2::DESeqDataSetFromTximport(txi.deseq, pdata, design = ~group)
  dds <- DESeq2::estimateSizeFactors(dds)
  vsd <- trans_fun(dds, blind = FALSE)

  return(vsd)
}


#' Construct expression set
#'
#' @param quants \code{DGEList} with RNA-seq counts.
#' @param fdata \code{data.table} returned from \code{\link{setup_fdata}}.
#' @param annot Character vector with ensembldb package name. e.g. \code{'EnsDb.Hsapiens.v94'}. Returned from \code{\link{get_ensdb_package}}.
#'
#' @return \code{ExpressionSet} with input arguments are in corresponding atributes.
#' @keywords internal
#' @export
#'
construct_eset <- function(quants, fdata, annot, txi.deseq = NULL) {
  # remove duplicate rows of counts
  rn <- row.names(quants$counts)

  txi.deseq <- txi.deseq[1:3]
  for (name in names(txi.deseq))
    colnames(txi.deseq[[name]]) <-
    paste(colnames(txi.deseq[[name]]), name, sep = '_')

  # workaround: R crashed from unique(mat) with GSE93624
  counts <- data.frame(quants$counts, rn, stringsAsFactors = FALSE, check.names = FALSE)

  mat <- data.table::data.table(counts,
                                txi.deseq$abundance,
                                txi.deseq$counts,
                                txi.deseq$length, key = 'rn')

  mat <- mat[!duplicated(counts), ]

  # merge exprs and fdata
  dt <- merge(fdata, mat, by.y = 'rn', by.x = data.table::key(fdata), all.y = TRUE, sort = FALSE)
  dt <- as.data.frame(dt)
  row.names(dt) <- make.unique(dt[[1]])

  # remove edgeR assigned group column
  pdata <- quants$samples
  pdata$group <- NULL

  # create environment with counts for limma and txi.deseq values for plots
  e <- new.env()
  e$exprs <- as.matrix(dt[, row.names(pdata), drop=FALSE])

  if (!is.null(txi.deseq)) {
    e$abundance <- as.matrix(dt[, paste0(row.names(pdata), '_abundance'), drop=FALSE])
    e$counts <- as.matrix(dt[, paste0(row.names(pdata), '_counts'), drop=FALSE])
    e$length <- as.matrix(dt[, paste0(row.names(pdata), '_length'), drop=FALSE])

    colnames(e$abundance) <- gsub('_abundance$', '', colnames(e$abundance))
    colnames(e$counts) <- gsub('_counts$', '', colnames(e$counts))
    colnames(e$length) <- gsub('_length$', '', colnames(e$length))
  }

  # seperate fdata and exprs and transfer to eset
  eset <- Biobase::ExpressionSet(e,
                                 phenoData=Biobase::AnnotatedDataFrame(pdata),
                                 featureData=Biobase::AnnotatedDataFrame(dt[, colnames(fdata), drop=FALSE]),
                                 annotation=annot)

  return(eset)
}

#' Setup feature annotation data
#'
#' @return \code{data.table} with columns \code{SYMBOL} and \code{ENTREZID_HS} corresponding to
#'   HGNC symbols and human entrez ids respectively.
#' @keywords internal
#' @export
#'
setup_fdata <- function(species = 'Homo sapiens', release = '94') {
  is.hs <- grepl('sapiens', species)
  if (species == 'Mus musculus') tx2gene <- tx2gene_mouse

  if (!grepl('sapiens|musculus', species)) {
    tx2gene <- get_tx2gene(species, release, columns = c("tx_id", "gene_name", "entrezid"))
  }

  # unlist entrezids
  fdata <- data.table::data.table(tx2gene)
  fdata <- fdata[, list(ENTREZID = as.character(unlist(entrezid)))
                 , by = c('tx_id', 'gene_name')]

  # add homologene
  fdata <- merge(unique(fdata), homologene, by = 'ENTREZID', all.x = TRUE, sort = FALSE)

  # add HGNC
  exist_homologues <- sum(!is.na(fdata$ENTREZID_HS)) != 0
  if (!exist_homologues) {
    fdata[, tx_id := NULL]
    fdata <- unique(fdata)

  } else {
    # where no homology, use original entrez id (useful if human platform):
    if (is.hs) {
      filt <- is.na(fdata$ENTREZID_HS)
      fdata[filt, "ENTREZID_HS"] <- fdata$ENTREZID[filt]
    }

    # map human entrez id to human symbol used by hs (for cmap/l1000 stuff)
    fdata$SYMBOL <- toupper(hs[fdata$ENTREZID_HS, SYMBOL_9606])
    fdata[, tx_id := NULL]

    if (is.hs) {
      # if gene_name (AnnotationHub) and SYMBOL (NCBI) differ use NCBI
      diff <- which(fdata$SYMBOL != fdata$gene_name)
      fdata$gene_name[diff] <- fdata$SYMBOL[diff]

      # if have gene_name but no SYMBOL use gene_name
      filt <- is.na(fdata$SYMBOL)
      fdata$SYMBOL[filt] <- fdata$gene_name[filt]
    }

    fdata <- unique(fdata)

    # keep duplicate gene_name -> SYMBOL and gene only if non NA
    discard <- fdata[, .I[any(!is.na(SYMBOL)) & is.na(SYMBOL)], by=gene_name]$V1
    if (length(discard)) fdata <- fdata[-discard]
  }

  data.table::setkeyv(fdata, 'gene_name')
  return(fdata)
}

#' Import kallisto quants
#'
#' @inheritParams setup_fdata
#' @inheritParams load_seq
#'
#' @return \code{DGEList} with length scaled counts.
#' @keywords internal
#' @export
#'
import_quants <- function(data_dir, species = 'Homo sapiens', release = '94') {

  if (!grepl('sapiens', species)) {
    tx2gene <- get_tx2gene(species, release, columns = c("tx_id", "gene_name", "entrezid"))
  }

  # don't ignoreTxVersion if dots in tx2gene
  ignore <- TRUE
  if (any(grepl('[.]', tx2gene$tx_id))) ignore <- FALSE

  # import quants using tximport
  # using limma::voom for differential expression (see tximport vignette)
  pkg_version <- get_pkg_version('kallisto')
  qdir <- paste('kallisto', pkg_version, 'quants', sep = '_')
  qdirs <- list.files(file.path(data_dir, qdir))
  quants_paths <- file.path(data_dir, qdir, qdirs, 'abundance.h5')

  # use folders as names (used as sample names)
  names(quants_paths) <- qdirs

  # import limma for differential expression analysis
  txi.limma <- tximport::tximport(quants_paths, tx2gene = tx2gene, type = 'kallisto',
                                  ignoreTxVersion = ignore, countsFromAbundance = 'lengthScaledTPM')

  quants <- edgeR::DGEList(txi.limma$counts)
  quants <- edgeR::calcNormFactors(quants)

  # import DESeq2 for exploratory analysis (plots)
  txi.deseq <- tximport::tximport(quants_paths, tx2gene = tx2gene, type = 'kallisto',
                                  ignoreTxVersion = ignore, countsFromAbundance = 'no')


  return(list(quants = quants, txi.deseq = txi.deseq))
}


#' Get ensembldb package name
#'
#' @inheritParams load_seq
#'
#' @keywords internal
#' @return Character vector with ensembldb package name. e.g. \code{'EnsDb.Hsapiens.v94'}.
#' @export
#'
get_ensdb_package <- function(species, release) {
  ensdb_species    <- strsplit(species, ' ')[[1]]
  ensdb_species[1] <- toupper(substr(ensdb_species[1], 1, 1))

  ensdb_package <- paste('EnsDb', paste0(ensdb_species, collapse = ''), paste0('v', release), sep='.')
  return(ensdb_package)
}

#' Get transcript to gene map.
#'
#' @inheritParams load_seq
#' @param columns Character vector of columns from ensdb package to return or \code{'list'} to print available options.
#'
#' @return \code{data.frame} with columns \code{tx_id}, \code{gene_name}, and \code{entrezid}
#' @export
#'
#' @keywords internal
#'
get_tx2gene <- function(species = 'Homo sapiens', release = '94', columns = c("tx_id", "gene_name", "entrezid", "gene_id", "seq_name", "description")) {
  is.installed('ensembldb', level = 'error')

  # load EnsDb package
  ensdb_package <- get_ensdb_package(species, release)
  if (!require(ensdb_package, character.only = TRUE)) {
    build_ensdb(species, release)
    require(ensdb_package, character.only = TRUE)
  }

  # for printing available columns in EnsDb package
  if (columns[1] == 'list') {
    print(ensembldb::listColumns(get(ensdb_package)))
    return(NULL)
  }

  # map from transcripts to genes
  tx2gene <- ensembldb::transcripts(get(ensdb_package), columns=columns, return.type='data.frame')
  tx2gene[tx2gene == ""] <- NA
  tx2gene <- tx2gene[!is.na(tx2gene$gene_name), ]

  if ('description' %in% colnames(tx2gene))
    tx2gene$description <- gsub(' \\[Source.+?\\]', '', tx2gene$description)

  return(tx2gene)
}


#' Build ensembldb annotation package.
#'
#' @inheritParams load_seq
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # build ensembldb annotation package for human
#' build_ensdb()
#'
build_ensdb <- function(species = 'Homo sapiens', release = '94') {
  # packages in suggests for faster install time
  suggests <- c('ensembldb', 'AnnotationHub', 'AnnotationDbi')
  is.installed(suggests, level = 'error')

  # store ensembl databases in built package
  ensdb_dir <- 'EnsDb'
  unlink('EnsDb', recursive = TRUE)
  dir.create(ensdb_dir)

  # format is genus_species in multiple other functions but not here
  species <- gsub('_', ' ', species)

  # generate new ensembl database from specified release
  ah <- AnnotationHub::AnnotationHub()
  ahDb <- AnnotationHub::query(ah, pattern = c(species, "EnsDb", release))

  if (!length(ahDb)) stop('Specified ensemble species/release not found in AnnotationHub.')

  ahEdb <- ahDb[[1]]

  ensembldb::makeEnsembldbPackage(AnnotationDbi::dbfile(ensembldb::dbconn(ahEdb)),
                                  '0.0.1', 'Alex Pickering <alexvpickering@gmail.com>',
                                  'Alex Pickering',
                                  ensdb_dir)

  # install new ensemble database
  ensdb_name <- list.files(ensdb_dir)
  ensdb_path <- file.path(ensdb_dir, ensdb_name)
  install.packages(ensdb_path, repos = NULL)

  # remove source files
  unlink(ensdb_dir, recursive = TRUE)
}
