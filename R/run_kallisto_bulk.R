#' Runs kallisto quantification for bulk samples.
#'
#' For pair-ended experiments, reads for each pair should be in a seperate file.
#'
#' @inheritParams get_kallisto_index
#' @param indices_dir Directory with kallisto indices. See \code{\link{build_kallisto_index}}.
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#' @param quant_meta Previous result of \code{\link{get_quant_meta}} (for fastqs downloaded by `GEOfastq`) or
#'   \code{\link{select_pairs}} (for non-public fastqs). Must contain column \code{'File Name'} and optionally
#'   \code{'Pair'} (rows with same value taken as paired), and \code{Replicate} (rows with same valued taken as replicates). \code{quant_meta} is used to bypass
#'   call to \code{select_pairs} GUI.
#' @param paired Boolean indicating if fastq files are paired. If \code{NULL} (default), fastqs are considered paired if they are non \code{NA}
#'   entries in \code{quant_meta$Pair}.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#' @param fl.mean Estimated average fragment length (only relevant for single-end reads). Default (\code{NULL}) uses 200.
#' @param fl.sd Estimated standard deviation of fragment length (only relevant for single-end reads). Default (\code{NULL}) uses 20.
#' @param updateProgress Used by drugseqr app to provide visual update of progress.
#'
#' @return NULL
#' @export
#'
#'
run_kallisto_bulk <- function(indices_dir, data_dir, quant_meta = NULL, paired = NULL, species = 'homo_sapiens', release = '94', fl.mean = NULL, fl.sd = NULL, updateProgress = NULL) {

  data_dir <- path.expand(data_dir)

  # default updateProgress and number of steps
  if (is.null(updateProgress)) updateProgress <- function(...) {NULL}

  # get index_path
  kallisto_version <- get_pkg_version('kallisto')
  index_path <- get_kallisto_index(indices_dir, species, release)

  if (is.null(quant_meta)) quant_meta <- select_pairs(data_dir)

  # save quants here
  quants_dir <- file.path(data_dir, paste('kallisto', kallisto_version, 'quants', sep = '_'))
  dir.create(quants_dir, showWarnings = FALSE)

  # if not supplied from app, then from select_pairs quant_meta
  if (is.null(paired)) paired <- sum(!is.na(quant_meta$Pair)) > 0

  # specific flags for single end experiments
  flags <- NULL
  if (!paired) {
    if (is.null(fl.mean)) {
      message('Single-end experiment but estimated average fragment length not provided. Setting to 200.')
      fl.mean <- 200
    }

    if (is.null(fl.sd)) {
      message('Single-end experiment but estimated standard deviation of fragment length not provided. Setting to 20.')
      fl.sd <- 20
    }
    flags <- c('--single', '-l', fl.mean, '-s', fl.sd)
  }

  # loop through fastq files and quantify
  fastq_files <- quant_meta$`File Name`
  if (is.null(fastq_files)) stop("File names should be in 'File Name' column")

  while (length(fastq_files)) {

    # grab next
    fastq_file <-fastq_files[1]
    row_num <- which(quant_meta$`File Name` == fastq_file)

    # include any pairs
    pair_num <- quant_meta$Pair[row_num]
    if (!is.null(pair_num) && !is.na(pair_num)) {
      fastq_file <- quant_meta[quant_meta$Pair %in% pair_num, 'File Name', drop = TRUE]
    }

    # include any replicates
    rep_num <- quant_meta$Replicate[row_num]
    if (!is.null(rep_num) && !is.na(rep_num)) {
      in.rep <- quant_meta$Replicate %in% rep_num
      fastq_file <- quant_meta[in.rep, 'File Name', drop = TRUE]

      # order by pairs
      pair_nums <- quant_meta$Pair[in.rep]
      if (!is.null(pair_nums) && !any(is.na(pair_nums))) {
        fastq_file <- fastq_file[order(pair_nums)]
      }
    }

    # update progress
    updateProgress(amount = length(fastq_file))

    # remove fastq_files for next quant loop
    fastq_files <- setdiff(fastq_files, fastq_file)

    # possibly spaces in file names
    fastq_path <- paste(shQuote(file.path(data_dir, fastq_file)), collapse = ' ')

    # save each sample in it's own folder
    # use the first file name in any replicates/pairs
    sample_name <- gsub('.fastq.gz$', '', fastq_file[1])
    out_dir <- file.path(quants_dir, sample_name)
    unlink(out_dir, recursive = TRUE)

    # run kallisto
    system2('kallisto',
            args=c('quant',
                   '-i', index_path,
                   '-o', shQuote(out_dir),
                   '-t', 6,
                   flags,
                   fastq_path))
  }
  return(NULL)
}


#' Get path to kallisto index
#'
#' @param indices_dir Path to folder with indices dir
#' @param species Character vector giving species
#' @param release Character vector of EnsDB release number. Release \code{'94'} is the default because HGNC symbols
#'   most closely overlap with those used to annotate CMAP02 and L1000.
#'
#' @return Path to kallisto index
#' @keywords internal
#'
get_kallisto_index <- function(indices_dir, species = 'homo_sapiens', release = '94') {
  kallisto_version <- get_pkg_version('kallisto')
  species <- gsub(' ', '_', tolower(species))

  index_path <- file.path(indices_dir, paste0('kallisto_', kallisto_version))
  index_path <- list.files(index_path, paste0(species, '.+?.cdna.all.release-', release, '_k31.idx'), full.names = TRUE)

  return(index_path)
}
