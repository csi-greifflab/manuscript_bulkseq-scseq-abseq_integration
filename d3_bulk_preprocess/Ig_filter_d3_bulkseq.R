library(tidyverse)
library(immunarch)

# Filter for only clonotypes with correct isotype, removing clonotypes with unidentified isotype and misidentified isotypes
path <- list.files("./1_Ig_all_unfiltered/", full.names = T)
file_names <- list.files("1_Ig_all_unfiltered/")
Ig_sorted_unfiltered <- lapply(path, function(x) read_tsv(x))
names(Ig_sorted_unfiltered) <- file_names

#Keep only clones that does not have an empty cdr3 & isotype
Ig_sorted_filtered <- lapply(Ig_sorted_unfiltered, function(x) x[which(is.na(x$aaSeqCDR3) == FALSE),])
Ig_sorted_filtered <- lapply(Ig_sorted_filtered, function(x) x[which(is.na(x$allCHitsWithScore) == FALSE),])

# Remove clones with light chain constant regions (dataset only contains IgH)
Ig_sorted_filtered_2 <- lapply(Ig_sorted_filtered, function(x) 
                          x[-str_which(x$allCHitsWithScore, "IGKC|IGLC"),])

# No light chain in some samples so add them back
Ig_sorted_filtered_2$BCP1_mem.txt <- Ig_sorted_filtered$BCP1_mem.txt
Ig_sorted_filtered_2$BCP1_naive.txt <- Ig_sorted_filtered$BCP1_naive.txt
Ig_sorted_filtered_2$BCP3_GC.txt <- Ig_sorted_filtered$BCP3_GC.txt
Ig_sorted_filtered_2$BCP3_mem.txt <- Ig_sorted_filtered$BCP3_mem.txt
Ig_sorted_filtered_2$BCP3_naive.txt <- Ig_sorted_filtered$BCP3_naive.txt
Ig_sorted_filtered_2$BCP4_mem.txt <- Ig_sorted_filtered$BCP4_mem.txt
Ig_sorted_filtered_2$BCP5_naive.txt <- Ig_sorted_filtered$BCP5_naive.txt
Ig_sorted_filtered_2$BCP6_all.txt <- Ig_sorted_filtered$BCP6_all.txt
Ig_sorted_filtered_2$BCP6_naive.txt <- Ig_sorted_filtered$BCP6_naive.txt
Ig_sorted_filtered_2$BCP8_all.txt <- Ig_sorted_filtered$BCP8_all.txt
Ig_sorted_filtered_2$BCP8_GC.txt <- Ig_sorted_filtered$BCP8_GC.txt
Ig_sorted_filtered_2$BCP8_naive.txt <- Ig_sorted_filtered$BCP8_naive.txt
Ig_sorted_filtered_2$BCP9_GC.txt <- Ig_sorted_filtered$BCP9_GC.txt
Ig_sorted_filtered_2$BCP9_mem.txt <- Ig_sorted_filtered$BCP9_mem.txt
Ig_sorted_filtered_2$BCP9_naive.txt <- Ig_sorted_filtered$BCP9_naive.txt


# Write out the files as tsv
map2(.x = Ig_sorted_filtered_2, .y = file_names, function(x,y)
     write_tsv(x, file = paste0("./2_Ig_all_filtered/", y)))

# Load with immunarch
data_d3_bulk <- repLoad("./2_Ig_all_filtered/")

# Recalculate the clone count, clonal proportions, and add C gene information

for (i in 1:length(data_d3_bulk$data)){
  data_d3_bulk$data[[i]]$Clones <- Ig_sorted_filtered_2[[i]]$uniqueUMICount
  data_d3_bulk$data[[i]]$Best.C <- str_split(str_split(data_d3_bulk$data[[i]]$C.name, ",", simplify = T)[,1], "\\*", simplify = T)[,1]
  data_d3_bulk$data[[i]]$Proportion <- data_d3_bulk$data[[i]]$Clones / sum(data_d3_bulk$data[[i]]$Clones)
}

# Save the pre-processed immunarch table as RDS file
saveRDS(data_d3_bulk, "./data_d3_bulk.RDS")
