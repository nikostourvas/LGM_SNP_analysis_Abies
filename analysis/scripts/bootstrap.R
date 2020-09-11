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

### Calculate parameter
# library(hierfstat)
# x <- poppr2hierfstat_out(obj_SSR, "Hexp")
# # OR
# x <- 1 / (1 - (basic.stats(obj_SSR)[["Hs"]]))
# # OR
# x <- allelic.richness(obj_SSR, min.n=44)
# x <- x[["Ar"]]
# 
# 

### bootstraping a df of a parameter
boot.param.sd <- function(x, nboot=1000){
        
        if(is.vector(x) == TRUE){
                boot.mean <- function(vec, i){
                        return(mean(vec[i]))
                }
                res <- boot(data = x,
                            statistic = boot.mean,
                            R = nboot)
                return(sd(res$t))
                
        } else {
                boot.mean.df <- function(df, i){
                        df2 <- df[i,]
                        return(colMeans(df2))
                }
                res <- boot(data = x,
                            statistic = boot.mean.df,
                            R=nboot)
                res <- apply(res$t, 2, sd)
                res <- data.frame(sd = res)
                rownames(res) <- colnames(x)
                return(res)
        }
}