setwd("~/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline")
library('sleuth')

#modify directory todata
base_dir <- "."

sample_id <- dir(file.path(base_dir, 'singlecellquant2cluster0'))
sample_id <- sample_id[1:40]
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "singlecellquant2cluster0", id))
covariate <- rep('one', length(sample_id))

sample_id2 <- dir(file.path(base_dir, 'singlecellquant2cluster1'))
sample_id2 <- sample_id2[1:40]
kal_dirs2 <- sapply(sample_id2, function(id) file.path(base_dir, "singlecellquant2cluster1", id))
covariate2 <- rep('two', length(sample_id2))

sample_id3 <- dir(file.path(base_dir, 'singlecellquant2cluster2'))
sample_id3 <- sample_id3[1:40]
kal_dirs3 <- sapply(sample_id3, function(id) file.path(base_dir, "singlecellquant2cluster2", id))
covariate3 <- rep('three', length(sample_id3))


sample_id <- append(sample_id, sample_id2)
sample_id <- append(sample_id, sample_id3)

kal_dirs <- append(kal_dirs, kal_dirs2)
kal_dirs <- append(kal_dirs, kal_dirs3)

covariate <-append(covariate, covariate2)
covariate <- append(covariate, covariate3)

fragments <- seq(from =1, to=1, length.out=length(sample_id))

#modify directory of design matrix
s2c <- data.frame(sample=sample_id, cluster = covariate, fragments = fragments)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(s2c, ~ cluster, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#modify sleuth_live call?
sleuth_live(so)
