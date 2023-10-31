library(readr)
# Filter for only clonotypes with correct isotype, removing clonotypes with unidentified isotype and misidentified isotypes
IgA_i5 <- read_delim("1_Ig_all_unfiltered/IgA_i5.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgA_i5_f <- IgA_i5[grep("^IGHA",IgA_i5$allCHitsWithScore),]
IgA_i5_f <- IgA_i5_f[which(is.na(IgA_i5_f$aaSeqCDR3) == FALSE),]
write_tsv(IgA_i5_f,"2_Ig_all_filtered/IgA_i5.txt")

IgD_i6 <- read_delim("1_Ig_all_unfiltered/IgD_i6.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgD_i6_f <- IgD_i6[grep("^IGHD",IgD_i6$allCHitsWithScore),]
IgD_i6_f <- IgD_i6_f[which(is.na(IgD_i6_f$aaSeqCDR3) == FALSE),]
write_tsv(IgD_i6_f,"2_Ig_all_filtered/IgD_i6.txt")

IgE_i7 <- read_delim("1_Ig_all_unfiltered/IgE_i7.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgE_i7_f <- IgE_i7[grep("^IGHE",IgE_i7$allCHitsWithScore),]
IgE_i7_f <- IgE_i7_f[which(is.na(IgE_i7_f$aaSeqCDR3) == FALSE),]
write_tsv(IgE_i7_f,"2_Ig_all_filtered/IgE_i7.txt")

IgG_i8 <- read_delim("1_Ig_all_unfiltered/IgG_i8.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgG_i8_f <- IgG_i8[grep("^IGHG",IgG_i8$allCHitsWithScore),]
IgG_i8_f <- IgG_i8_f[which(is.na(IgG_i8_f$aaSeqCDR3) == FALSE),]
write_tsv(IgG_i8_f,"2_Ig_all_filtered/IgG_i8.txt")

IgK_i10 <- read_delim("1_Ig_all_unfiltered/IgK_i10.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgK_i10_f <- IgK_i10[grep("^IGKC",IgK_i10$allCHitsWithScore),]
IgK_i10_f <- IgK_i10_f[which(is.na(IgK_i10_f$aaSeqCDR3) == FALSE),]
write_tsv(IgK_i10_f,"2_Ig_all_filtered/IgK_i10.txt")

IgL_i11 <- read_delim("1_Ig_all_unfiltered/IgL_i11.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgL_i11_f <- IgL_i11[grep("^IGLC",IgL_i11$allCHitsWithScore),]
IgL_i11_f  <- IgL_i11_f[which(is.na(IgL_i11_f$aaSeqCDR3) == FALSE),]
write_tsv(IgL_i11_f,"2_Ig_all_filtered/IgL_i11.txt")

IgM_i9 <- read_delim("1_Ig_all_unfiltered/IgM_i9.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
IgM_i9_f <- IgM_i9[grep("^IGHM",IgM_i9$allCHitsWithScore),]
IgM_i9_f  <- IgM_i9_f[which(is.na(IgM_i9_f$aaSeqCDR3) == FALSE),]
write_tsv(IgM_i9_f,"2_Ig_all_filtered/IgM_i9.txt")