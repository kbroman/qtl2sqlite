#' Write genotype probabilities to database
#'
#' Write genotype probabilities to database
#'
#' @param db Either a filename for a SQLite database, or a connection to such a database.
#' @param probs Genotype probabilities, as calculated by \code{\link[qtl2geno]{calc_genoprob}}.
#' @param map Optional map of marker positions (a list of vectors of positions).
#' @param quiet If FALSE, print some information about progress.
#'
#' @return None.
#'
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable
#' @importFrom RSQLite SQLite dbGetQuery
#' @export
write_probs <-
    function(db, probs, map=NULL, quiet=FALSE)
{
    # check inputs
    if(!is.null(map)) match_probs_map(probs, map)
    nind <- vapply(probs, nrow, 1)
    if(length(unique(nind)) != 1)
        stop("probs has different numbers of individuals on different chromosomes")
    ind <- lapply(probs, rownames)
    if(length(ind) > 1) { # more than one chromosome
        for(i in seq(along=ind)[-1]) {
            if(any(ind[[1]] != ind[[i]]))
                stop("probs has different row names on different chromosomes")
        }
    }

    if(is.character(db)) { # filename, so connect (and disconnect on exit)
        if(file.exists(db)) {
            unlink(db) # if file exists, remove it
            if(!quiet) message("Writing over file ", db)
            warning("Writing over file ", db)
        }

        db <- DBI::dbConnect(RSQLite::SQLite(), db)
        on.exit(DBI::dbDisconnect(db))
    }

    # write chromosome table
    is_x_chr <- attr(probs, "is_x_chr")
    if(is.null(is_x_chr)) {
        warning("Missing is_x_chr attribute; assuming all are autosomes")
        is_x_chr <- rep(FALSE, length(probs))
    }
    chr_tab <- data.frame(chr=names(probs),
                          is_x_chr=is_x_chr,
                          chr_index=seq(along=probs),
                          stringsAsFactors=FALSE)
    DBI::dbWriteTable(db, "chr", chr_tab, row.names=FALSE,
                      overwrite=TRUE)

    # write alleles attribute
    alleles <- attr(probs, "alleles")
    if(!is.null(alleles)) {
        allele_tab <- data.frame(alleles=alleles,
                                 allele_index=seq(along=alleles),
                                 stringsAsFactors=FALSE)
        DBI::dbWriteTable(db, "alleles", allele_tab, row.names=FALSE,
                          overwrite=TRUE)
    }

    # write attribute table
    crosstype <- attr(probs, "crosstype")
    if(is.null(crosstype)) crosstype <- ""
    alleleprobs <- attr(probs, "alleleprobs")
    if(is.null(alleleprobs)) alleleprobs <- FALSE
    attr_tab <- data.frame(crosstype=crosstype,
                           alleleprobs=alleleprobs)
    DBI::dbWriteTable(db, "attributes", attr_tab, row.names=FALSE,
                      overwrite=TRUE)


    # write markers table with map information, if provided
    if(!is.null(map)) {
        map <- map_list_to_df(map)
        if(length(unique(map$marker)) != nrow(map))
            stop("Marker names are not unique")

        DBI::dbWriteTable(db, "markers", map, row.names=FALSE, overwrite=TRUE)
        dbGetQuery(db, "CREATE INDEX markers_pos ON markers (chr, pos)")
    }
    else { # no map, so just make table with chr and marker
        mn <- lapply(probs, function(a) dimnames(a)[[3]])
        nmar <- vapply(probs, function(a) dim(a)[3], 1)
        map <- data.frame(chr=rep(names(probs), nmar),
                          marker=unlist(mn),
                          marker_index=unlist(lapply(mn, function(a) seq(along=a))),
                          stringsAsFactors=FALSE)

        DBI::dbWriteTable(db, "markers", map, row.names=FALSE, overwrite=TRUE)
    }
    dbGetQuery(db, "CREATE INDEX markers_chr ON markers (chr)")

    # write individual table
    ind_tab <- data.frame(ind=rownames(probs[[1]]),
                          ind_index=1:nrow(probs[[1]]))
    DBI::dbWriteTable(db, "ind", ind_tab, row.names=FALSE,
                      overwrite=TRUE)

    for(i in seq(along=probs)) {
        # write geno table
        chr <- names(probs)[i]
        if(!quiet) message("Writing probs for chr ", chr)
        geno <- colnames(probs[[i]])
        n_geno <- length(geno)
        geno_tab <- data.frame(chr=rep(chr, n_geno),
                               geno=geno,
                               geno_index=seq(along=geno),
                               stringsAsFactors=FALSE)
        DBI::dbWriteTable(db, "geno", geno_tab, row.names=FALSE,
                          append=TRUE)

        probs_tab <- probs_array2tab(probs, chr=i)
        print(dimnames(probs[[i]])[[3]])
        DBI::dbWriteTable(db, "probs", probs_tab, row.names=FALSE,
                          append=TRUE)
    }
    dbGetQuery(db, "CREATE INDEX geno_chr ON geno (chr)")
    dbGetQuery(db, "CREATE INDEX probs_marker ON probs (marker)")

    invisible(NULL)
}


# check that probs and map conform
match_probs_map <-
    function(probs, map, map_name="map")
{
    nmar_probs <- vapply(probs, function(a) dim(a)[[3]], 1)

    if(!all(nmar_probs == vapply(map, length, 1)))
        stop("probs and ", map_name, " have different numbers of markers")
    if(!all(unlist(lapply(probs, function(a) dimnames(a)[[3]])) ==
            unlist(lapply(map, names))))
        stop("probs and ", map_name, " have different marker names")

    TRUE
}
