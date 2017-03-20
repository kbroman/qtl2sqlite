# convert map from list to data frame
#   include_index: if TRUE, include "marker_index" column with numeric indices within chr
#   row.names: if TRUE, use marker names as row names
#   if marker_column is NULL, don't include a "marker" column
#       (and include markers as row names even if row.names=FALSE)
map_list_to_df <-
    function(map_list, chr_column="chr", pos_column="pos", marker_column="marker",
             include_index=TRUE, row.names=FALSE)
{
    nmar <- vapply(map_list, length, 1) # no. markers per chromosome

    markers <- unlist(lapply(map_list, names))

    result <- data.frame(chr=rep(names(map_list), nmar),
                         pos=unlist(map_list),
                         marker=markers,
                         stringsAsFactors=FALSE)
    if(include_index)
        result$marker_index <- unlist(lapply(map_list, function(a) seq(along=a)))
    if(row.names)
        rownames(result) <- markers

    names(result)[1] <- chr_column
    names(result)[2] <- pos_column
    if(is.null(marker_column) && row.names)
        result <- result[,-3,drop=FALSE]
    else if(!is.null(marker_column))
        names(result)[3] <- marker_column

    result
}

# convert map from df to list
map_df_to_list <-
    function(map, chr_column="chr", pos_column="cM", marker_column="marker",
             Xchr=c("x", "X"))
{
    if(is.null(marker_column)) {
        marker_column <- "qtl2tmp_marker"
        map[,marker_column] <- rownames(marker_column)
    }
    if(!(marker_column %in% colnames(map)))
        stop('Column "', marker_column, '" not found.')
    if(!(chr_column %in% colnames(map)))
        stop('Column "', chr_column, '" not found.')
    if(!(pos_column %in% colnames(map)))
        stop('Column "', pos_column, '" not found.')

    marker <- map[,marker_column]

    chr <- map[,chr_column]
    uchr <- unique(chr)
    pos <- map[,pos_column]

    result <- split(pos, factor(chr, levels=uchr))
    marker <- split(marker, factor(chr, levels=uchr))
    for(i in seq(along=result))
        names(result[[i]]) <- marker[[i]]

    is_x_chr <- rep(FALSE, length(result))
    names(is_x_chr) <- names(result)
    if(!is.null(Xchr)) {
        Xchr_used <- Xchr %in% names(is_x_chr)
        if(any(Xchr_used)) {
            Xchr <- Xchr[Xchr_used]
            is_x_chr[Xchr] <- TRUE
        }
    }
    attr(result, "is_x_chr") <- is_x_chr

    result
}
