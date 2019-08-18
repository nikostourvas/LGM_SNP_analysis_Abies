
snp <- read.csv(file = "../data/Genotyping-1841.025-03 Grid_reformated.csv", 
                header = T, 
                na.strings = c("?", "Uncallable", "Bad", "Missing")
                # ,stringsAsFactors = T
                , check.names = F # default is TRUE, and changes the dashes to dots - which created problems
)
# snp




library(tidyverse)
snp <- snp %>% 
        mutate(Genotype = str_replace_all(Genotype, "^AB", "GR_AB_A_")) %>% #GR_Adult
        mutate(Genotype = str_replace_all(Genotype, "^RAB", "GR_AB")) %>% #GR_Regen
        
        mutate(Genotype = str_replace_all(Genotype, "_A_", "_A_1_")) %>% 
        
        mutate(Genotype = str_replace_all(Genotype, "NR", "NR_"))


### Create a df following the guidelines of the loci format

snp_loci_format <- as.data.frame(sapply(snp, gsub, pattern=':', replacement='/') )
# snp_loci_format <- as.data.frame(sapply(snp_loci_format, gsub, pattern=".", replacement="_"))
snp_loci_format <- snp_loci_format[,-1]
rownames(snp_loci_format) <- snp[,1]

# snp_loci_format


# convert NAs to -/-
dartR_df <- as.matrix(snp_loci_format)
dartR_df[is.na(dartR_df)] <- "-/-"
dartR_df <- as.data.frame(dartR_df)

write.csv(dartR_df, "../data/dartR_input.csv")

