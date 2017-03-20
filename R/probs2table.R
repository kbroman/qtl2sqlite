# convert probs for a chromosome to a big table
probs_array2tab <-
    function(probs, chr=1)
{
    d <- dim(probs[[chr]])
    dn <- dimnames(probs[[chr]])
    data.frame(ind=rep(dn[[1]], d[2]*d[3]),
               geno=rep(rep(dn[[2]], each=d[1]), d[3]),
               marker=rep(dn[[3]], each=d[1]*d[2]),
               prob=as.numeric(probs[[chr]]),
               stringsAsFactors=FALSE)

}

# convert probs (for one chromosome) from a table back to a 3-d array
# *** critical that the rows are in the right order! (marker than genotype than ind'l)
probs_tab2array <-
    function(table, ind=NULL, geno=NULL, marker=NULL)
{
    if(is.null(ind))
        ind <- unique(table$ind)
    if(is.null(geno))
        geno <- unique(table$geno)
    if(is.null(marker))
        marker <- unique(table$marker)

    p <- array(table$prob, dim=c(length(ind), length(geno), length(marker)))
    dimnames(p) <- list(ind, geno, marker)

    p
}
