#' Transform GEO metadata for quantification
#'
#' @param srp_meta data.frame of GEO metadata from
#'   \code{\link[GEOfastq]{crawl_gsms}}
#' @inheritParams run_kallisto_bulk
#'
#' @return data.frame needed for parameter \code{quant_meta} in
#'   \code{\link{run_kallisto_bulk}}
#' @export
#' @examples
#'
#' # example srp_meta
#' srp_meta <- data.frame(
#'     run = "SRR1",
#'     gsm_name = "GSM4875733",
#'     library_layout = "PAIRED"
#' )
#'
#' # example paired fastq files
#' data_dir <- tempdir()
#' file.create(file.path(data_dir, c("SRR1_1.fastq.gz", "SRR1_2.fastq.gz")))
#' quant_meta <- get_quant_meta(srp_meta, data_dir)
get_quant_meta <- function(srp_meta, data_dir) {
    fastq_files <- list.files(
        data_dir, paste0(srp_meta$run, "(_[12])?.fastq.gz$", collapse = "|")
    )

    srr_names <- gsub(".fastq.gz$", "", fastq_files)
    srr_names <- gsub("_[12]$", "", srr_names)

    # check for paired fastq files
    layout <- unique(srp_meta$library_layout)
    if (length(layout) != 1) stop("Multiple library layouts not implemented")

    pair <- NA
    paired <- layout == "PAIRED"

    if (paired) {
        pair <- as.numeric(as.factor(srr_names))
        has.pair <- duplicated(pair) | duplicated(pair, fromLast = TRUE)
        stopifnot(all(has.pair))
    }

    # check for GSMs split between multiple fastq files
    rep <- NA
    reped <- any(duplicated(srp_meta$gsm_name))

    if (reped) {
        rep <- as.numeric(as.factor(srp_meta[srr_names, "gsm_name"]))
        has.rep <- duplicated(rep) | duplicated(rep, fromLast = TRUE)
        rep[!has.rep] <- NA
    }


    data.frame(
        `File Name` = fastq_files,
        Pair = pair,
        Replicate = rep,
        run = srr_names, stringsAsFactors = FALSE, check.names = FALSE
    )
}
