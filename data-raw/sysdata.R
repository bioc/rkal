tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')
tx2gene_mouse <- readRDS('data-raw/tx2gene/tx2gene_mouse.rds')
homologene <- readRDS('data-raw/homologene/homologene.rds')
hs <- readRDS('data-raw/homologene/hs.rds')


usethis::use_data(tx2gene, tx2gene_mouse, homologene, internal = TRUE, overwrite = TRUE)
usethis::use_data(hs, internal = FALSE, overwrite = TRUE)
