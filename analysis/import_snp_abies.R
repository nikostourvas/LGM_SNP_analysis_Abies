# customized function
import.snp.abies <- function(threshold_loci, threshold_ind, maf){

snp <- read.csv("../data/Genotyping-1841.025-03 Grid_reformated.csv", 
                header = T, 
                na.strings = c("?", "Uncallable", "Bad", "Missing")
                # ,stringsAsFactors = T
                , check.names = F # default is TRUE, and changes the dashes to dots - which created problems
)
# snp



 ## Transform sample names 
         #### This will allow hierarchical analysis when applicable 
         #### Country / Species / Pop / Plot 
        
library(tidyverse)

snp <- snp %>% 
        mutate(Genotype = str_replace_all(Genotype, "^AB", "GR_AB_A_")) %>% #GR_Adult
        mutate(Genotype = str_replace_all(Genotype, "^RAB", "GR_AB")) %>% #GR_Regen
        
        mutate(Genotype = str_replace_all(Genotype, "_A_", "_A_1_")) %>% 
        
        mutate(Genotype = str_replace_all(Genotype, "NR", "NR_"))


 ## Remove samples that consistently do not amplify 
        
# 
# snp_fil <- snp %>% 
#         filter(snpa != "Bad")
# 
# snp_fil


 ### Create a df following the guidelines of the loci format 
        
snp_loci_format <- as.data.frame(sapply(snp, gsub, pattern=':', replacement='/') )
# snp_loci_format <- as.data.frame(sapply(snp_loci_format, gsub, pattern=".", replacement="_"))
snp_loci_format <- snp_loci_format[,-1]
rownames(snp_loci_format) <- snp[,1]

# snp_loci_format



## Add a second column for population

pop1 <- replicate(124, "DE_AAL_A_1")
pop2 <- replicate(16, "DE_AAL_NR_1")
pop3 <- replicate(16, "DE_AAL_NR_2")
pop4 <- replicate(16, "DE_AAL_NR_3")
pop5 <- replicate(16, "DE_AAL_NR_4")
pop6 <- replicate(124, "GR_AB_A_1")
pop7 <- replicate(20, "GR_AB_NR_1")
pop8 <- replicate(18, "GR_AB_NR_2")
pop9 <- replicate(16, "GR_AB_NR_3")
pop10 <- replicate(10, "GR_AB_NR_4")
pop11 <- replicate(124, "SI_AAL_A_1")
pop12 <- replicate(16, "SI_AAL_NR_1")
pop13 <- replicate(16, "SI_AAL_NR_2")
pop14 <- replicate(16, "SI_AAL_NR_3")
pop15 <- replicate(16, "SI_AAL_NR_4")

pop <- c(pop1, pop2, pop3, pop4, pop5, pop6,
         pop7, pop8, pop9, pop10, pop11,
         pop12, pop13, pop14, pop15)

# pop

snp_loci_format <- add_column(snp_loci_format, pop, .before = "contig00241-160")

# snp_loci_format


 ## Create genind object 
        
library(pegas)

data <- as.loci(snp_loci_format, 
                col.pop = 1
                ,allele.sep = "/")
# data

obj_origin <- loci2genind(data,)
# obj_origin


 ### stratify data set 
        
strata_df <- as.data.frame(snp_loci_format$pop)
colnames(strata_df) <- "strata"
strata_df <- separate(strata_df, col = strata, sep="_", 
                      into = c("Country", 
                               "Species", 
                               "Pop", 
                               "Plot"))
strata(obj_origin) <- strata_df

setPop(obj_origin) <- ~Country/Pop
# obj_origin


 # Data filtering 
        
         ### Remove monomorphic loci 
        
library(poppr)
obj <- informloci(obj_origin, MAF = 0)


 ### Filter out missing data 

threshold_loci <- threshold_loci
threshold_ind <- threshold_ind

obj <- missingno(obj, type = "loci", cutoff = threshold_loci, quiet = T)

obj <- missingno(obj, type = "genotypes", cutoff = threshold_ind, quiet = T)


 ### Remove uninformative loci 
        
library(poppr)

maf <- maf

obj <- informloci(obj, MAF = maf)




return(obj)
}




# by pop ------------------------------------------------------------------

# customized function
# the only difference with import.snp.abies is the handling
# of loci missing data
# Here it is performed by pop

import.snp.abies2 <- function(threshold_loci, threshold_ind, maf){
  
  snp <- read.csv("../data/Genotyping-1841.025-03 Grid_reformated.csv", 
                  header = T, 
                  na.strings = c("?", "Uncallable", "Bad")
                  # ,stringsAsFactors = T
                  , check.names = F # default is TRUE, and changes the dashes to dots - which created problems
  )
  # snp
  
  
  
  ## Transform sample names 
  #### This will allow hierarchical analysis when applicable 
  #### Country / Species / Pop / Plot 
  
  library(tidyverse)
  
  snp <- snp %>% 
    mutate(Genotype = str_replace_all(Genotype, "^AB", "GR_AB_A_")) %>% #GR_Adult
    mutate(Genotype = str_replace_all(Genotype, "^RAB", "GR_AB")) %>% #GR_Regen
    
    mutate(Genotype = str_replace_all(Genotype, "_A_", "_A_1_")) %>% 
    
    mutate(Genotype = str_replace_all(Genotype, "NR", "NR_"))
  
  
  ## Remove samples that consistently do not amplify 
  
  # 
  # snp_fil <- snp %>% 
  #         filter(snpa != "Bad")
  # 
  # snp_fil
  
  
  ### Create a df following the guidelines of the loci format 
  
  snp_loci_format <- as.data.frame(sapply(snp, gsub, pattern=':', replacement='/') )
  # snp_loci_format <- as.data.frame(sapply(snp_loci_format, gsub, pattern=".", replacement="_"))
  snp_loci_format <- snp_loci_format[,-1]
  rownames(snp_loci_format) <- snp[,1]
  
  # snp_loci_format
  
  
  
  ## Add a second column for population
  
  pop1 <- replicate(124, "DE_AAL_A_1")
  pop2 <- replicate(16, "DE_AAL_NR_1")
  pop3 <- replicate(16, "DE_AAL_NR_2")
  pop4 <- replicate(16, "DE_AAL_NR_3")
  pop5 <- replicate(16, "DE_AAL_NR_4")
  pop6 <- replicate(124, "GR_AB_A_1")
  pop7 <- replicate(20, "GR_AB_NR_1")
  pop8 <- replicate(18, "GR_AB_NR_2")
  pop9 <- replicate(16, "GR_AB_NR_3")
  pop10 <- replicate(10, "GR_AB_NR_4")
  pop11 <- replicate(124, "SI_AAL_A_1")
  pop12 <- replicate(16, "SI_AAL_NR_1")
  pop13 <- replicate(16, "SI_AAL_NR_2")
  pop14 <- replicate(16, "SI_AAL_NR_3")
  pop15 <- replicate(16, "SI_AAL_NR_4")
  
  pop <- c(pop1, pop2, pop3, pop4, pop5, pop6,
           pop7, pop8, pop9, pop10, pop11,
           pop12, pop13, pop14, pop15)
  
  # pop
  
  snp_loci_format <- add_column(snp_loci_format, pop, .before = "contig00241-160")
  
  # snp_loci_format
  
  
  ## Create genind object 
  
  library(pegas)
  
  data <- as.loci(snp_loci_format, 
                  col.pop = 1
                  ,allele.sep = "/")
  # data
  
  obj_origin <- loci2genind(data,)
  # obj_origin
  
  
  ### stratify data set 
  
  strata_df <- as.data.frame(snp_loci_format$pop)
  colnames(strata_df) <- "strata"
  strata_df <- separate(strata_df, col = strata, sep="_", 
                        into = c("Country", 
                                 "Species", 
                                 "Pop", 
                                 "Plot"))
  strata(obj_origin) <- strata_df
  
  setPop(obj_origin) <- ~Country/Pop
  # obj_origin
  
  
  # Data filtering 
  
  ### Remove monomorphic loci 
  
  library(poppr)
  obj <- informloci(obj_origin, MAF = 0)
  
  
  ### Filter out missing data 
  
  threshold_loci <- threshold_loci
  threshold_ind <- threshold_ind
  
  obj <- missingno(obj, type = "loci", cutoff = threshold_loci, quiet = T)
  
  obj <- missingno(obj, type = "genotypes", cutoff = threshold_ind, quiet = T)
  
  
  
  info <- info_table(obj, type = "missing", plot = F, plotlab = F)
  
  # Remove loci with specified % of missing data by pop
  miss <- as.data.frame(t(info))
  miss$loc <- rownames(miss)
  
  miss_res <- miss %>% 
    filter_if(is.double, any_vars(. > threshold_loci))
  
  miss_res$loc
  
  # create list of loci to keep
  keeploc <- setdiff(locNames(obj), miss_res$loc)
  
  # filter loci in genind object
  obj <- obj[loc = keeploc]
  
  
  
  
  ### Remove uninformative loci 
  
  library(poppr)
  
  maf <- maf
  
  obj <- informloci(obj, MAF = maf)
  
  
  
  
  return(obj)
}



