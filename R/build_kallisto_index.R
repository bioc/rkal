#' Download ensembl transcriptome and build index for kalliso quantification
#'
#' This index is used for RNA seq quantification.
#'
#' @param indices_dir Path to directory to create indices folder in for storing kallisto indices.
#' @param species The species. Default is \code{homo_sapiens}.
#' @param release ensembl release. Default is \code{94} (latest in release for AnnotationHub -
#'   needs to match with \code{\link{build_ensdb}}).
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build kallist index for humans
#' indices_dir <- '~/kallisto'
#' build_kallisto_index(indices_dir)
#'
build_kallisto_index <- function(indices_dir, species = 'homo_sapiens', release = '94') {

  kallisto_version <- get_pkg_version('kallisto')
  indices_dir <- file.path(indices_dir, paste0('kallisto_', kallisto_version))

  if (!dir.exists(indices_dir))
    dir.create(indices_dir, recursive = TRUE)

  # construct ensembl url for transcriptome
  ensembl_species <- gsub(' ', '_', tolower(species))
  ensembl_release <- paste0('release-', release)
  ensembl_url <- paste0('ftp://ftp.ensembl.org/pub/', ensembl_release, '/fasta/', ensembl_species, '/cdna/')

  # get list of all files
  handle <- curl::new_handle(dirlistonly=TRUE)
  con <- curl::curl(ensembl_url, "r", handle)
  tbl <- utils::read.table(con, stringsAsFactors=FALSE)
  close(con)

  # get transcripts cdna.all file
  ensembl_fasta <- grep('cdna.all.fa.gz$', tbl[[1]], value = TRUE)
  ensembl_url <- paste0(ensembl_url, ensembl_fasta)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(ensembl_url, ensembl_fasta)

  # build index
  index_fname <- gsub('fa.gz$', paste0(ensembl_release, '_k31.idx'), ensembl_fasta)
  index_fname <- tolower(index_fname)
  tryCatch(system2('kallisto', args=c('index',
                                      '-i', index_fname,
                                      ensembl_fasta)),
           error = function(err) {err$message <- 'Is kallisto installed and on the PATH?'; stop(err)})

  unlink(ensembl_fasta)
  setwd(work_dir)
}

#' Get version of salmon/kallisto from system command.
#'
#' @param type Either \code{'salmon'} or \code{'kallisto'}.
#'
#' @return Version of salmon/kallisto.
#' @keywords internal
#'
#' @examples
get_pkg_version <- function(type) {
  # possibly use older salmon with version appended to executable name
  if (type == 'salmon') {
    version <- system(paste(type, '--version'), intern = TRUE)
    version <- gsub('^salmon ', '', version)

  } else if (type == 'kallisto') {
    version <- system(paste(type, 'version'), intern = TRUE)
    version <- gsub('^kallisto, version ', '', version)

  } else {
    stop('type must be either salmon or kallisto')
  }

  return(version)
}