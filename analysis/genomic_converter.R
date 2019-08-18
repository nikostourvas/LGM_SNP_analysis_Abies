library(radiator)

test <- genomic_converter(obj, output = "genlight")
glPlot(test$genlight, posi="topleft")

temp <- seppop(test$genlight) 
myFreq <- lapply(temp, glMean)
myFreq <- lapply(myFreq, function(x){
        c(x, 1-x)
})


par(mfrow = c(3, 2)) 

hist(myFreq[[1]], proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="DE_A", nclass=20)
hist(myFreq[[2]], proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="DE_NR", nclass=20)
hist(myFreq[[3]], proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="GR_A", nclass=20)
hist(myFreq[[4]], proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="GR_NR", nclass=20)
hist(myFreq[[5]], proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="SI_A", nclass=20)
hist(myFreq[[6]], proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="SI_NR", nclass=20)

temp <- seppop(test$genlight) 
myFreq <- lapply(temp, glMean)
hist(myFreq[[3]], proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)



#-----------------------------------------
# pcadapt
obj_pcadapt <- genomic_converter(obj, output = "pcadapt")

library("pcadapt")
library("qvalue")

pca_genotype <- read.pcadapt(obj_pcadapt$pcadapt$genotype.matrix)
K <- 25



x1 <- pcadapt(pca_genotype, K = K) # if it fails to run, copy this
# line and run it in the console
p1 <- plot(x1, option = "screeplot") # 4 groups seems to be the correct value


loci <- genind2loci(obj)
library(pegas)
p2 <- plot(x1, option = "scores", pop = loci[,1])

setPop(obj) <- ~Country/Pop      


K <- 4

x <- pcadapt(pca_genotype, K = K, min.maf = 0)

summary(x) # numerical quantities obtained after performing a PCA



p3 <- plot(x, option = "manhattan")
