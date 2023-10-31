library(tidyverse)
library(seqinr)

# Read the processed MiXCR sequences
path_bulk <- list.files("../../../20210921_bulkseq/Analysis/20210921_bulkseq_mixcr_v4100/2_Ig_all_filtered/", full.names = T)
data_bulk <- map(path_bulk, function(x) read_tsv(x))
# exclude the bulk IgE samples
data_bulk <- data_bulk[-3]
names(data_bulk) <- c("Sample_0_IgA_bulk", "Sample_0_IgD_bulk", "Sample_0_IgG_bulk", "Sample_0_IgK_bulk", "Sample_0_IgL_bulk", "Sample_0_IgM_bulk")
data_bulk_hc <- data_bulk[c(1:3,6)]
data_bulk_lc <- data_bulk[c(4,5)]

path_sc <- list.files("../../../20210714_scseq/Analysis/20210714_scseq_mixcr_v4100_merged/6_Ig_all_clonotypes_isotypes_single/", full.names = T)
names_sc <- str_remove_all(list.files("../../../20210714_scseq/Analysis/20210714_scseq_mixcr_v4100_merged/6_Ig_all_clonotypes_isotypes_single/", full.names = F), ".txt")
data_sc <- map(path_sc, function(x) read_tsv(x))
names(data_sc) <- names_sc
data_sc_hc <- data_sc[-str_which(names(data_sc), "IgK|IgL")]
data_sc_lc <- data_sc[str_which(names(data_sc), "IgK|IgL")]

vdj_bulk_hc <- map(data_bulk_hc, function(x) as.list(toupper(str_remove_all(x$aaSeqImputedVDJRegion, "_"))))
vdj_bulk_lc <- map(data_bulk_lc, function(x) as.list(toupper(str_remove_all(x$aaSeqImputedVDJRegion, "_"))))
vdj_sc_hc <- map(data_sc_hc, function(x) as.list(toupper(str_remove_all(x$aaSeqImputedVDJRegion, "_"))))
vdj_sc_lc <- map(data_sc_lc, function(x) as.list(toupper(str_remove_all(x$aaSeqImputedVDJRegion, "_"))))

clone_id_bulk_hc <-  map(data_bulk_hc, function(x) x$cloneId)
clone_id_bulk_lc <-  map(data_bulk_lc, function(x) x$cloneId)
clone_id_sc_hc <-  map(data_sc_hc, function(x) x$cloneId)
clone_id_sc_lc <-  map(data_sc_lc, function(x) x$cloneId)

# add the clone ID number to the names of the vdj sequence
for (i in 1:length(vdj_bulk_hc)){
  names(vdj_bulk_hc[[i]]) <- clone_id_bulk_hc[[i]]
}

for (i in 1:length(vdj_bulk_lc)){
  names(vdj_bulk_lc[[i]]) <- clone_id_bulk_lc[[i]]
}

for (i in 1:length(vdj_sc_hc)){
  names(vdj_sc_hc[[i]]) <- clone_id_sc_hc[[i]]
}

for (i in 1:length(vdj_sc_lc)){
  names(vdj_sc_lc[[i]]) <- clone_id_sc_lc[[i]]
}

vdj2_bulk_hc <- vdj_bulk_hc %>% list_flatten(name_spec = "{outer}_{inner}")
vdj2_bulk_lc <- vdj_bulk_lc %>% list_flatten(name_spec = "{outer}_{inner}")
vdj2_sc_hc <- vdj_sc_hc %>% list_flatten(name_spec = "{outer}_{inner}")
vdj2_sc_lc <- vdj_sc_lc %>% list_flatten(name_spec = "{outer}_{inner}")
#vdj2 <- c(vdj2_bulk, vdj2_sc)

# Identify duplicates VDJ and rename them to include all clonotypes ID with that VDJ
# When running MiXCR, don't seperate clones by C name
# use distinct() on a df
# identify shared clones between bulk and sc as well

# write the vdj sequences to fasta
write.fasta(vdj2_bulk_hc, names = names(vdj2_bulk_hc), file.out = "./fasta/d1_vdj_bulk_hc.fasta", as.string = TRUE)
write.fasta(vdj2_bulk_lc, names = names(vdj2_bulk_lc), file.out = "./fasta/d1_vdj_bulk_lc.fasta", as.string = TRUE)
write.fasta(vdj2_sc_hc, names = names(vdj2_sc_hc), file.out = "./fasta/d1_vdj_sc_hc.fasta", as.string = TRUE)
write.fasta(vdj2_sc_lc, names = names(vdj2_sc_lc), file.out = "./fasta/d1_vdj_sc_lc.fasta", as.string = TRUE)
# extract cdr3 sequence
cdr3_bulk <- map(data_bulk, function(x) x$aaSeqCDR3)
cdr3_sc <- map(data_sc, function(x) x$aaSeqCDR3)
cdr3 <- c(cdr3_bulk, cdr3_sc)

# write the vdj sequences to a tsv file
df_ref <- tibble(Clonotype = names(vdj2),
                 Sample = paste(str_split_i(names(vdj2), "_", 1), str_split_i(names(vdj2), "_", 2), str_split_i(names(vdj2), "_", 3), str_split_i(names(vdj2), "_", 4), sep = "_"),
                 Clone_ID = str_split(names(vdj2), "_", simplify = T)[,5],
                 VDJ_seq = unlist(vdj2),
                 CDR3_seq = unlist(cdr3),
                 Isotype = str_split(Sample, "_", simplify = T)[,3],
                 Method = str_split(Sample, "_", simplify = T)[,4])

write_tsv(df_ref, "./df_ref.tsv")
