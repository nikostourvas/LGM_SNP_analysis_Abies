priv_all <- function(obj){

priv <- poppr::private_alleles(obj, report = "data.frame")

if(is.data.frame(priv)){
        
        p1<-ggplot(priv) + geom_tile(aes(x = population,
                                         y = allele,
                                         fill = count)) +
                ggtitle("Private alleles per population") +
                scale_fill_viridis_c()
        
}

setPop(obj) <- ~Country

priv <- poppr::private_alleles(obj, 
                               report = "data.frame")


if(is.data.frame(priv)){
        
        p2<-ggplot(priv) + geom_tile(aes(x = population,
                                         y = allele,
                                         fill = count)) +
                ggtitle("Private alleles per cohort") +
                scale_fill_viridis_c()
        
}

setPop(obj) <- ~Country/Pop


        
setPop(obj) <- ~Country

priv_count <- poppr::private_alleles(obj, 
                                     report = "data.frame",
                                     count.alleles = F)
priv_count$maf <- maf

# ggplot(priv_count, aes(x=population, y=count)) +
#   geom_bar(stat="identity")

setPop(obj) <- ~Country/Pop



list <- list(p1, p2, priv_count)

}