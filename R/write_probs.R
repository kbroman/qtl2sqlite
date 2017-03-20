#' Write genotype probabilities to database
#'
#' Write genotype probabilities to database
#'
#' @param db Either a filename for a SQLite database, or a connection to such a database.
#' @param probs Genotype probabilities, as calculated by \code{\link[qtl2geno]{calc_genoprob}}.
#'
#' @return None.
#'
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable
#' @importFrom RSQLite SQLite dbGetQuery
#' @export
write_probs <-
    function(db, probs, pmap=NULL, gmap=NULL)
{

    if(is.character(db)) { # filename, so connect (and disconnect on exit)
        if(file.exists(db)) {
            unlink(db) # if file exists, remove it
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
                      overwrite=TRUE, append=FALSE)

    # write alleles attribute
    alleles <- attr(probs, "alleles")
    if(!is.null(alleles)) {
        allele_tab <- data.frame(alleles=alleles,
                                 allele_index=seq(along=alleles),
                                 stringsAsFactors=FALSE)
        DBI::dbWriteTable(db, "alleles", allele_tab, row.names=FALSE,
                          overwrite=TRUE, append=FALSE)
    }

    # write attribute table
    crosstype <- attr(probs, "crosstype")
    if(is.null(crosstype)) crosstype <- ""
    alleleprobs <- attr(probs, "alleleprobs")
    if(is.null(alleleprobs)) alleleprobs <- FALSE
    attr_tab <- data.frame(crosstype=crosstype,
                           alleleprobs=alleleprobs)
    DBI::dbWriteTable(db, "attributes", attr_tab, row.names=FALSE,
                      overwrite=TRUE, append=FALSE)


    # write map table, if provided
    if(!is.null(pmap)) { # physical map included
        match_probs_map(probs, pmap, "pmap")

        map <- map_list_to_df(pmap, pos_column="pos_Mbp")
        if(length(unique(map$marker)) != nrow(map))
            stop("Marker names are not unique")

        if(!is.null(gmap)) { # genetic map also included
            gmap_df <- map_list_to_df(gmap, pos="pos_cM")
            if(nrow(map) != nrow(gmap_df))
                stop("gmap and pmap have different numbers of markers")
            if(!all(gmap$marker == pmap$marker) ||
               !all(gmap$chr == pmap$chr) ||
               !all(gmap$marker_index == pmap$marker_index))
                stop("gmap and pmap have different marker/chr names")
            map$pos_cM <- gmap_df$pos_cM
        }
    }
    else if(!is.null(gmap)) { # just the genetic map
        match_probs_map(probs, gmap, "gmap")

        map <- map_list_to_df(pmap, pos_column="pos_cM")
        if(length(unique(map$marker)) != nrow(map))
            stop("Marker names are not unique")
    }
    if(!is.null(gmap) && !is.null(pmap)) {
        DBI::dbWriteTable(db, "map", map, row.names=FALSE, overwrite=TRUE, append=FALSE)
        if(!is.null(pmap))
            dbGetQuery(db, "CREATE INDEX pos_Mbp ON map (chr, pos_Mbp)")
        if(!is.null(gmap))
            dbGetQuery(db, "CREATE INDEX pos_cM ON map (chr, pos_cM)")
    }

    for(i in seq(along=probs)) {
    }

    # dbWriteTable(db, "table", tab)
    # dbGetQuery(db, "CREATE INDEX index_name ON table (var1, var2)")
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
