ic <- function(obj, P){
        na_by_locus <- poppr2hierfstat_out(obj, "allele")
        na <- table_out(obj, na_by_locus, "na")
        na <- na[nrow(na),1]
        # Ng
        ng <- na * (na+1) / 2
        # IC
        ic <- nLoc(obj) * P * ng
return(ic)
}