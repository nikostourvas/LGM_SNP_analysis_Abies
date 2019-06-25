barplot.N.alleles <- function(obj) {
div <- summary(obj)

all <- as.data.frame(div$pop.n.all)
all <- all %>% 
        mutate(population = rownames(all)) %>% 
        rename(count = "div$pop.n.all")

all$count <- as.integer(all$count)

# str(all)

ggplot(all, aes(x=population, y=count) ) +
        geom_bar(stat = "identity", 
                 color = "black", fill = "white",
                 width = .8) +
        geom_text(aes(label=count, vjust = 1.6)) +
        ggtitle("Total Number of alleles per population")

}