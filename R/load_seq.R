.datatable.aware <- TRUE

#' Load bulk RNA-Seq data into an ExpressionSet.
#'
#' @param data_dir Directory with raw and quantified RNA-Seq files.
#' @param species Character vector indicating species. Genus and species should
#'   be space separated, not underscore. Default is \code{Homo sapiens}.
#' @param release EnsemblDB release. Should be same as used in
#'   \code{\link{build_kallisto_index}}.
#' @param load_saved If TRUE (default) and a saved \code{ExpressionSet} exists,
#'   will load from disk.
#' @param save_eset If TRUE (default) and either \code{load_saved} is
#'   \code{FALSE} or a saved \code{ExpressionSet} does not exist,
#'   then an ExpressionSet will be saved to disk. Will overwrite if
#'   already exists.
#'
#' @return \code{ExpressionSet} with attributes/accessors:
#' \itemize{
#'   \item \code{sampleNames} from names of raw RNA-Seq files (excluding
#'     .fastq.gz suffix).
#'   \item \code{annotation} Character vector of annotation package used.
#'   \item \code{exprs} Length scaled counts generated from abundances for use
#'     in \code{\link[limma]{voom}} (see
#'     \code{vignette("tximport", package = "tximport")}).
#'   \item \code{abundance, counts, length} accessed e.g. by
#'     \code{assayDataElement(eset, 'length')}.
#'   Imported for exploratory data analysis with DESeq2 variance stabilization
#'     transforms.
#'   \item \code{phenoData} added columns:
#'     \itemize{
#'       \item \code{lib.size} library size from
#'         \code{\link[edgeR]{calcNormFactors}}.
#'       \item \code{norm.factors} library normalization factors from
#'         \code{\link[edgeR]{calcNormFactors}}.
#'     }
#' }
#'
#' @export
#' @seealso \link{get_vsd}
#' @examples
#'
#' library(tximportData)
#' example <- system.file("extdata", "kallisto_boot",
#'     "ERR188021",
#'     package = "tximportData"
#' )
#'
#' # setup to mirror expected folder structure
#' data_dir <- tempdir()
#' qdir <- paste("kallisto", get_pkg_version(), "quants", sep = "_")
#' qdir <- file.path(data_dir, qdir)
#' unlink(qdir, recursive = TRUE)
#' dir.create(qdir)
#' file.copy(example, qdir, recursive = TRUE)
#'
#' # construct and annotate eset
#' eset <- load_seq(data_dir)
load_seq <- function(data_dir, species = "Homo sapiens", release = "94",
    load_saved = TRUE, save_eset = TRUE) {

    # check if already have
    eset_path <- file.path(data_dir, "eset.rds")
    if (load_saved & file.exists(eset_path)) {
        return(readRDS(eset_path))
    }

    # import quants and filter low counts
    q <- import_quants(data_dir, species = species, release = release)

    # construct eset
    annot <- drugseqr.data::get_ensdb_package(species, release)
    fdata <- setup_fdata(species, release)
    eset <- construct_eset(q$quants, fdata, annot, q$txi.deseq)

    # save eset and return
    if (save_eset) saveRDS(eset, eset_path)
    return(eset)
}

#' Load ExpressionSet from ARCHS4 h5 file
#'
#' ARCHS4 h5 files are available at https://maayanlab.cloud/archs4/download.html
#'
#' @param archs4_file Path to human_matrix.h5 file.
#' @param gsm_names Character vector of GSM names of samples to load.
#' @param eset_path Path to load saved eset from or to save eset to.
#' @param load_saved Should a previously saved eset be loaded? Requires
#'   \code{eset_path} to be specified.
#' @inheritParams load_seq
#'
#'
#' @return ExpressionSet
#' @export
#'
#' @examples
#' archs4_file <- "/path/to/human_matrix_v*.h5"
#' gsm_names <- c("GSM3190508", "GSM3190509", "GSM3190510", "GSM3190511")
#' # eset <- load_archs4_seq(archs4_file, gsm_names)
load_archs4_seq <- function(archs4_file, gsm_names, species = "Homo sapiens",
    release = "94", load_saved = TRUE, eset_path = NULL) {
    saved_eset <- !is.null(eset_path) && file.exists(eset_path)
    if (saved_eset & load_saved) {
        return(readRDS(eset_path))
    }

    samples <- as.character(rhdf5::h5read(
        archs4_file,
        "meta/samples/geo_accession"
    ))
    genes <- as.character(rhdf5::h5read(archs4_file, "meta/genes/genes"))

    # use ARCHS4
    sample_locations <- which(samples %in% gsm_names)
    if (length(sample_locations) < length(gsm_names)) {
        warning(
            "Only ", length(sample_locations), " of ", length(gsm_names),
            " GSMs present."
        )
    }

    # extract gene expression from compressed data
    counts <- t(rhdf5::h5read(archs4_file, "data/expression",
        index = list(sample_locations, seq_along(genes))
    ))
    rhdf5::H5close()

    counts[counts < 0] <- 0
    rownames(counts) <- genes
    colnames(counts) <- samples[sample_locations]

    quants <- edgeR::calcNormFactors(edgeR::DGEList(counts))
    fdata <- setup_fdata(species, release)
    annot <- drugseqr.data::get_ensdb_package(species, release)
    eset <- construct_eset(quants, fdata, annot)

    if (!is.null(eset_path)) saveRDS(eset, eset_path)
    return(eset)
}

#' Get variance stabilized data for exploratory data analysis
#'
#' @param eset ExpressionSet loaded with \link{load_seq}.
#'   Requires group column in \code{pData(eset)} specifying sample groupings.
#' @param rlog_cutoff Sample number above which will use
#'   \code{\link[DESeq2]{vst}} instead of \code{\link[DESeq2]{rlog}}.
#'   Default is 50.
#'
#' @return \code{DESeqTransform} with variance stabilized expression data.
#' @export
#' @examples
#'
#' # generate example
#' y <- matrix(rnbinom(10000, mu = 5, size = 2), ncol = 4)
#' row.names(y) <- paste0("gene", 1:2500)
#' quants <- edgeR::DGEList(counts = y)
#'
#' fdata <- data.table::data.table(gene_name = row.names(y), key = "gene_name")
#' annot <- drugseqr.data::get_ensdb_package("Homo sapiens", "94")
#'
#' eset <- construct_eset(quants, fdata, annot)
#' eset$group <- factor(c("t", "t", "c", "c"))
#' vsd <- get_vsd(eset)
get_vsd <- function(eset, rlog_cutoff = 50) {
    trans_fun <- if (ncol(eset) > rlog_cutoff) DESeq2::vst else DESeq2::rlog
    els <- Biobase::assayDataElementNames(eset)
    pdata <- Biobase::pData(eset)

    if (all(c("abundance", "counts", "length") %in% els)) {
        txi.deseq <- list(
            countsFromAbundance = "no",
            abundance = Biobase::assayDataElement(eset, "abundance"),
            counts = Biobase::assayDataElement(eset, "counts"),
            length = Biobase::assayDataElement(eset, "length")
        )

        dds <- DESeq2::DESeqDataSetFromTximport(txi.deseq,
            pdata,
            design = ~group
        )
    } else {
        # this is e.g. for eset from load_archs4_seq
        dds <- DESeq2::DESeqDataSetFromMatrix(Biobase::exprs(eset),
            pdata,
            design = ~group
        )
    }
    dds <- DESeq2::estimateSizeFactors(dds)
    vsd <- trans_fun(dds, blind = FALSE)

    return(vsd)
}


#' Construct expression set
#'
#' @param quants \code{DGEList} with RNA-seq counts.
#' @param fdata \code{data.table} returned from \code{\link{setup_fdata}}.
#' @param annot Character vector with ensembldb package name. e.g.
#'   \code{'EnsDb.Hsapiens.v94'}. Returned from
#'   \code{\link[drugseqr.data]{get_ensdb_package}}.
#' @param txi.deseq Optional \code{DGElist} returned by \code{import_quants}. If
#'   specified, assays 'counts', 'abundance', and 'length' will be present in
#'   the returned \code{ExpressionSet}.
#'
#' @return \code{ExpressionSet} with \code{quants$counts} stored in the 'exprs'
#'  assay, \code{quants$samples} stored in the sample metadata, and \code{fdata}
#'  stored in the feature data.
#' @export
#' @examples
#'
#' # generate example
#' y <- matrix(rnbinom(10000, mu = 5, size = 2), ncol = 4)
#' row.names(y) <- paste0("gene", 1:2500)
#' quants <- edgeR::DGEList(counts = y)
#'
#' fdata <- data.table::data.table(gene_name = row.names(y), key = "gene_name")
#' annot <- drugseqr.data::get_ensdb_package("Homo sapiens", "94")
#'
#' eset <- construct_eset(quants, fdata, annot)
construct_eset <- function(quants, fdata, annot, txi.deseq = NULL) {
    # remove duplicate rows of counts
    rn <- row.names(quants$counts)

    txi.deseq <- txi.deseq[c("abundance", "counts", "length")]
    for (name in names(txi.deseq)) {
        colnames(txi.deseq[[name]]) <-
            paste(colnames(txi.deseq[[name]]), name, sep = "_")
    }

    # workaround: R crashed from unique(mat) with GSE93624
    counts <- data.frame(quants$counts, rn,
        stringsAsFactors = FALSE, check.names = FALSE
    )

    mat <- data.table::data.table(counts,
        txi.deseq$abundance,
        txi.deseq$counts,
        txi.deseq$length,
        key = "rn"
    )

    mat <- mat[!duplicated(counts), ]

    # merge exprs and fdata
    dt <- merge(fdata, mat,
        by.y = "rn", by.x = data.table::key(fdata),
        all.y = TRUE, sort = FALSE
    )
    dt <- as.data.frame(dt)
    row.names(dt) <- make.unique(dt[[1]])

    # remove edgeR assigned group column
    pdata <- quants$samples
    pdata$group <- NULL

    # create environment with counts for limma and txi.deseq values for plots
    e <- new.env()
    e$exprs <- as.matrix(dt[, row.names(pdata), drop = FALSE])

    if (!is.null(txi.deseq)) {
        e$abundance <- as.matrix(
            dt[, paste0(row.names(pdata), "_abundance"), drop = FALSE]
        )
        e$counts <- as.matrix(
            dt[, paste0(row.names(pdata), "_counts"), drop = FALSE]
        )
        e$length <- as.matrix(
            dt[, paste0(row.names(pdata), "_length"), drop = FALSE]
        )

        colnames(e$abundance) <- gsub("_abundance$", "", colnames(e$abundance))
        colnames(e$counts) <- gsub("_counts$", "", colnames(e$counts))
        colnames(e$length) <- gsub("_length$", "", colnames(e$length))
    }

    # seperate fdata and exprs and transfer to eset
    eset <- Biobase::ExpressionSet(
        e,
        phenoData = Biobase::AnnotatedDataFrame(pdata),
        featureData = Biobase::AnnotatedDataFrame(dt[, colnames(fdata),
            drop = FALSE
        ]),
        annotation = annot
    )

    return(eset)
}

#' Setup feature annotation data
#'
#' @return \code{data.table} with columns \code{SYMBOL} and \code{ENTREZID_HS}
#'   corresponding to HGNC symbols and human entrez ids respectively.
#' @importFrom rlang :=
#' @inheritParams build_kallisto_index
#' @export
#' @examples
#' setup_fdata()
setup_fdata <- function(species = "Homo sapiens", release = "94") {
    # for R CMD check
    entrezid <- tx_id <- SYMBOL_9606 <- NULL
    SYMBOL <- ENTREZID <- gene_name <- .I <- NULL

    is.hs <- grepl("sapiens", species)
    tx2gene <- drugseqr.data::load_tx2gene(species, release)

    # unlist entrezids
    fdata <- data.table::data.table(tx2gene)
    fdata <- fdata[, list(ENTREZID = as.character(unlist(entrezid))),
        by = c("tx_id", "gene_name")
    ]

    # add homologene
    homologene <- readRDS(system.file("extdata",
        "homologene.rds",
        package = "drugseqr.data"
    ))

    fdata <- merge(unique(fdata), homologene,
        by = "ENTREZID",
        all.x = TRUE, sort = FALSE
    )

    # add HGNC
    exist_homologues <- sum(!is.na(fdata$ENTREZID_HS)) != 0
    if (!exist_homologues) {
        fdata[, tx_id := NULL]
        fdata <- unique(fdata)
    } else {
        hs <- readRDS(system.file("extdata",
            "hs.rds",
            package = "drugseqr.data"
        ))

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

        # keep duplicate SYMBOL and ENTREZID only if non-NA
        discard <- fdata[,
            .I[(any(!is.na(SYMBOL)) & is.na(SYMBOL)) |
                (any(!is.na(ENTREZID)) & is.na(ENTREZID))],
            by = gene_name
        ]$V1

        if (length(discard)) fdata <- fdata[-discard]
    }

    data.table::setkeyv(fdata, "gene_name")
    return(fdata)
}

#' Import kallisto quants
#'
#' @inheritParams setup_fdata
#' @inheritParams load_seq
#'
#' @return list with items \itemize{
#'   \item{quants}{\code{DGEList} from \code{tximport} with argument
#'     \code{countsFromAbundance = 'lengthScaledTPM'}}
#'   \item{txi.deseq}{\code{DGEList}from \code{tximport} with argument
#'     \code{countsFromAbundance = 'no'}}
#' }
#' @keywords internal
#' @export
#' @examples
#' library(tximportData)
#' example <- system.file("extdata", "kallisto_boot",
#'     "ERR188021",
#'     package = "tximportData"
#' )
#'
#' # setup to mirror expected folder structure
#' data_dir <- tempdir()
#' qdir <- paste("kallisto", get_pkg_version(), "quants", sep = "_")
#' qdir <- file.path(data_dir, qdir)
#' unlink(qdir, recursive = TRUE)
#' dir.create(qdir)
#' file.copy(example, qdir, recursive = TRUE)
#' quants <- import_quants(data_dir)
import_quants <- function(data_dir, species = "Homo sapiens", release = "94") {
    tx2gene <- drugseqr.data::load_tx2gene(species, release)

    # don't ignoreTxVersion if dots in tx2gene
    ignore <- TRUE
    if (any(grepl("[.]", tx2gene$tx_id))) ignore <- FALSE

    # import quants using tximport
    # using limma voom for differential expression (see tximport vignette)
    pkg_version <- get_pkg_version()
    qdir <- paste("kallisto", pkg_version, "quants", sep = "_")
    qdirs <- list.files(file.path(data_dir, qdir))
    quants_paths <- file.path(data_dir, qdir, qdirs, "abundance.h5")

    # use folders as names (used as sample names)
    names(quants_paths) <- qdirs

    # import limma for differential expression analysis
    txi.limma <- tximport::tximport(
        quants_paths,
        tx2gene = tx2gene, type = "kallisto",
        ignoreTxVersion = ignore, countsFromAbundance = "lengthScaledTPM"
    )

    quants <- edgeR::DGEList(txi.limma$counts)
    quants <- edgeR::calcNormFactors(quants)

    # import DESeq2 for exploratory analysis (plots)
    txi.deseq <- tximport::tximport(
        quants_paths,
        tx2gene = tx2gene, type = "kallisto",
        ignoreTxVersion = ignore, countsFromAbundance = "no"
    )


    return(list(quants = quants, txi.deseq = txi.deseq))
}
