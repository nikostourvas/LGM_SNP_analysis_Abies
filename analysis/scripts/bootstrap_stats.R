boot.over.loci <- function(obj, nboot=1000){

## Sep object by locus
sep_loci <- seploc(obj)

library(pegas)
sep_loci_pegas <- lapply(sep_loci, as.loci)
sep_loci_pegas <- lapply(sep_loci_pegas, function(x){x[, 2]})


## Create loci obj

obj_loci <- as.loci(obj)

obj_loci_base <- as.data.frame(obj_loci[, 1])
colnames(obj_loci_base) <- "population"


## Loops
samples <- list()
dfs <- list()
nboot <- nboot
for(i in 1:nboot){
        samples[[i]] <- sample(sep_loci_pegas, size=nLoc(obj), replace = TRUE)
}

# list to data.frame
dfs <- lapply(samples, function(x){do.call(cbind.data.frame, x)})

for(i in 1:nboot){
        dfs[[i]] <- cbind(obj_loci_base, dfs[[i]])
        dfs[[i]] <- as.loci(dfs[[i]])
}

## Create genind object
genind_objects <- lapply(dfs, loci2genind)

return(genind_objects)
}



###### basic statistics reported from poppr and their SEs
poppr2hierfstat_out <- function(obj, variable){
        
        obj_list <- seppop(obj)
        
        stats_poppr <- list()
        for(i in 1: length(obj_list)){
                stats_poppr[[i]] <- locus_table(obj_list[[i]], 
                                                information = F)
        }
        
        table_out <- list()
        for(i in 1:length(obj_list))
                table_out[[i]] <- stats_poppr[[i]][-nrow(stats_poppr[[1]]), 
                                                   variable]
        
        table_out <- as.matrix(as.data.frame(table_out))
        colnames(table_out) <- popNames(obj)
        
        return(table_out)
}




boot.ci <- function(obj, sim_data, statistic, min.ar=NA){
        
if (statistic == "Hexp"){        
        
        uHe_boot_by_loc <- lapply(sim_data, 
                                      poppr2hierfstat_out, "Hexp")
        
        boot <- lapply(uHe_boot_by_loc, colMeans, na.rm=T)
        boot <- do.call(rbind.data.frame, boot)
        
        # uhe_boot_sd_A <- sd(uHe_boot[,1])
        # uhe_boot_sd_NR <- sd(uHe_boot[,2])
                
                
        # empirical data set value        
        empirical <- poppr2hierfstat_out(obj, "Hexp")
        empirical <- colMeans(empirical)


} else if (statistic == "ne"){
        library(hierfstat)
        ne_by_locus_Hs <- lapply(sim_data, function(x) {
                1 / (1 - (basic.stats(x)[["Hs"]])) 
        }
        )
        
        boot <- lapply(ne_by_locus_Hs, colMeans, na.rm=T)
        boot <- do.call(rbind.data.frame, boot)
        
        # ne_sd_A <- sd(ne_boot[,1])
        # ne_sd_NR <- sd(ne_boot[,2])
        
        # empirical data set value
        empirical <- 1 / (1 - (basic.stats(obj)[["Hs"]]))
        empirical <- colMeans(empirical)
        

} else if (statistic == "Ar"){
        library(hierfstat)
        ar_by_locus <- lapply(sim_data, allelic.richness, min.n=min.ar)
        ar_by_locus <- lapply(ar_by_locus, function(x){x[["Ar"]]})
        
        boot <- lapply(ar_by_locus, colMeans, na.rm=T)
        boot <- do.call(rbind.data.frame, boot)
        
        # ar_sd_A <- sd(ar_boot[,1])
        # ar_sd_NR <- sd(ar_boot[,2])
        
        # empirical data set value
        empirical <- allelic.richness(obj, min.n = min.ar)
        empirical <- colMeans(empirical$Ar)

}
        
        ci <- list()
        for(i in 1:nPop(obj)){
                deltastar = boot[,i] - empirical[i]
                d = quantile(deltastar, c(0.025, 0.975))
                ci[[i]] <- empirical[i] - c(d[2], d[1])
        }
        
        ci <- do.call(rbind.data.frame, ci)
        colnames(ci) <- c("lower", "upper")
        ci$pop <- popNames(obj)
        
return(ci)
}