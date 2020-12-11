test_that("get_quant_meta pairs paired GEO fastq.gz files", {

  srp_meta <- data.frame(run = 'SRR1',
                        gsm_name = 'GSM4875733',
                        library_layout = 'PAIRED')

  # example paired fastq files
  data_dir <- tempdir()
  file.create(file.path(data_dir, c('SRR1_1.fastq.gz', 'SRR1_2.fastq.gz')))
  quant_meta <- get_quant_meta(srp_meta, data_dir)

  expect_equal(quant_meta$Pair, c(1,1))
})

test_that("get_quant_meta doesn't pair single-end GEO fastq.gz files", {

  srp_meta <- data.frame(run = 'SRR1',
                         gsm_name = 'GSM4875733',
                         library_layout = 'SINGLE')

  # example paired fastq files
  data_dir <- tempdir()
  file.create(file.path(data_dir, c('SRR1_1.fastq.gz', 'SRR1_2.fastq.gz')))
  quant_meta <- get_quant_meta(srp_meta, data_dir)

  expect_equal(quant_meta$Pair, c(NA,NA))
})
