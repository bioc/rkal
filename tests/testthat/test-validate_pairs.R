test_that("get_fastq_id1s works", {
    fastq_paths <- system.file("testdata", "SRR12960930.400.fastq.gz", package = "rkal")
    fastq_id1s <- get_fastq_id1s(fastq_paths)

    expect_equal(unname(fastq_id1s), "@SRR12960930.1 1/1")
})


test_that("detect_paired works for single-end fastq.gz file", {
    fastq_paths <- system.file("testdata", "SRR12960930.400.fastq.gz", package = "rkal")
    fastq_id1s <- get_fastq_id1s(fastq_paths)


    expect_false(detect_paired(fastq_id1s))
})

test_that("detect_paired works for paired-end fastq.gz file", {
    fastq_dir <- system.file("testdata", package = "rkal")
    fastq_paths <- list.files(fastq_dir, "^SRR13202502.400_", full.names = TRUE)
    fastq_id1s <- get_fastq_id1s(fastq_paths)


    expect_true(detect_paired(fastq_id1s))
})


test_that("validate_pairs requires at least two samples to be paired", {
    pairs <- rep(NA, 4)
    rows <- c(1)
    reps <- rep(NA, 4)
    msg <- validate_pairs(pairs, rows, reps)
    expect_true(grepl("at least two", msg))
})


test_that("validate_reps requires at least two samples to be replicates", {
    pairs <- rep(NA, 4)
    rows <- c(1)
    reps <- rep(NA, 4)
    msg <- validate_reps(pairs, rows, reps)
    expect_true(grepl("at least two", msg))
})
