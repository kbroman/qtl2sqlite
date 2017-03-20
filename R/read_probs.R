#' Read genotype probabilities from database
#'
#' Read genotype probabilities from database
#'
#' @param db Either a filename for a SQLite database, or a connection to such a database.
#' @param chr Single chromosome to read.
#' @param pos Optional length-2 vector specifying an interval of positions to read.
#' @param markers Optional vector of marker names to read. If
#'     \code{pos} is provided, \code{marker} is ignored. The markers
#'     must all be on the same chromosome.
#'
#' @return Genotype probabilities, as an object of class \code{"calc_genoprob"}.
#'
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RSQLite SQLite dbGetQuery
#' @export
read_probs <-
    function(db, chr=NULL, pos=NULL, markers=NULL)
{
    if(is.character(db)) { # filename, so connect (and disconnect on exit)
        db <- DBI::dbConnect(RSQLite::SQLite(), db)
        on.exit(DBI::dbDisconnect(db))
    }

    if(is.null(chr) && is.null(markers))
        stop("Must specify either chr or markers")

    # read chromosome table
    chr_tab <- RSQLite::dbGetQuery(db, "SELECT chr, is_x_chr FROM chr ORDER BY chr_index")
    if(!is.null(chr)) {
        if(length(chr) != 1) {
            chr <- chr[1]
            warning('Argument "chr" must be a single character string.')
        }

        if(!(chr %in% chr_tab$chr))
            stop("Chromosome ", chr, " not in database")
    }
    is_x_chr <- as.logical(chr_tab$is_x_chr)
    names(is_x_chr) <- chr_tab$chr

    # attributes
    alleles <- RSQLite::dbGetQuery(db, "SELECT alleles FROM alleles")[,1]
    attrib <- RSQLite::dbGetQuery(db, "SELECT * FROM attributes")
    crosstype <- attrib$crosstype
    alleleprobs <- as.logical(attrib$alleleprobs)

    if(!is.null(pos)) {
        if(!is.null(markers))
            warning('Argument "markers" ignored if "pos" is provided')

        if(length(pos) != 2 || pos[1] > pos[2])
            stop('Argument "pos" should have length 2 with pos[1] <= pos[2]')

        if(!("pos" %in% dbListFields(db, "markers")))
            stop('markers table does not include position information, so argument "pos" can not be used.')

        markers <- RSQLite::dbGetQuery(db,
                              paste0('SELECT marker FROM markers WHERE chr=="', chr, '"',
                                     ' AND pos >= ', pos[1], ' AND pos <= ', pos[2],
                                     ' ORDER BY marker_index'))$marker

        if(length(markers) == 0)
            stop("No markers found in that interval")
    }
    else if(is.null(markers)) {
        markers <- RSQLite::dbGetQuery(db,
                             paste0('SELECT marker FROM markers WHERE chr=="', chr, '"',
                                    ' ORDER BY marker_index'))$marker
        if(length(markers) == 0)
            stop("No markers found on that chromosome")
    }

    # read ind
    ind <- RSQLite::dbGetQuery(db, "SELECT ind FROM ind ORDER BY ind_index")$ind

    # read probs
    pr <- RSQLite::dbGetQuery(db, paste0("SELECT probs.mat_index, probs.prob, probs.marker, markers.chr ",
                                "FROM probs, markers ",
                                "WHERE probs.marker == markers.marker ",
                                "AND markers.marker IN (",
                                paste0("'", markers, "'", collapse=","), ") ",
                                "ORDER BY markers.marker_index, probs.mat_index"))

    chr <- unique(pr$chr)
    if(length(chr) > 1)
        stop("markers on multiple chromosomes")

    # read geno for this chr
    geno <- RSQLite::dbGetQuery(db, paste0("SELECT geno FROM geno WHERE chr == '", chr, "'"))$geno

    # turn into array and add attributes
    pr <- list("1"=probs_tab2array(pr[,1:2,drop=FALSE], ind, geno, markers))
    names(pr) <- chr
    attr(pr, "crosstype") <- crosstype
    attr(pr, "is_x_chr") <- is_x_chr[chr]
    attr(pr, "alleles") <- alleles
    attr(pr, "alleleprobs") <- alleleprobs
    class(pr) <- c("calc_genoprob", "list")

    pr

}
