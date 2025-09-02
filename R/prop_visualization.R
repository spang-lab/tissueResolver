#' plots relative cell type fraction ratios between two groups
#'
#' @param props dataframe containing cell type proportions (output of cell_proportions)
#' @param grouping a dataframe containing "bulk_id" and "group" columns for group association of samples
#' @param groupA group that is in the numerator of the ratio
#' @param groupB group that is in the denominator of the ratio
#'
#' @return a ggplot object
#' 
#' @export
prop_visualization <- function(props, grouping, groupA, groupB) {
    if (!groupA %in% (grouping %>% dplyr::pull(group) %>% unique()) ||
        !groupB %in% (grouping %>% dplyr::pull(group) %>% unique())) {
        stop("given groups are not levels in the grouping file.")
    }
    if( ! "group" %in% names(grouping) ) {
      stop("no group in grouping df")
    }
    if( ! "bulk_id" %in% names(grouping) ) {
      stop("no bulk_id in grouping df")
    }


    groupvec <- grouping %>% dplyr::pull(group)
    names(groupvec) <- grouping %>% dplyr::pull(bulk_id)
    if( "weight_boot" %in% names(props) ){
        propdf <- props %>% tidyr::drop_na() %>%
            dplyr::select(-weight_boot)
    } else {
        propdf <- props %>% tidyr::drop_na()
    }
    if( ! "group" %in% names(props)) {
        propdf <- propdf %>% dplyr::inner_join(grouping, by="bulk_id")
    }
    propdf <- propdf %>%
        tidyr::pivot_wider(names_from = bulk_id, values_from=weight, values_fill = 0.0) %>%
        tidyr::pivot_longer(c(-celltype, -group), values_to = "weight", names_to="bulk_id") %>%
        dplyr::filter(groupvec[bulk_id] == group) %>% 
        dplyr::group_by(group, celltype) %>%
        dplyr::summarise(mweight = mean(weight), se=stats::sd(weight)/sqrt(length(weight))) %>%
        tidyr::pivot_wider(names_from = group, values_from=c(mweight, se)) %>%
        dplyr::mutate(
            frac = .data[[paste0("mweight_", groupA)]] / .data[[paste0("mweight_", groupB)]],
            se_frac = .data[[paste0("se_", groupA)]] / .data[[paste0("mweight_", groupA)]] +
                 .data[[paste0("se_", groupB)]] / .data[[paste0("mweight_", groupB)]]
            )
    return(
        ggplot2::ggplot(propdf, ggplot2::aes(x=celltype, y=frac, ymin=frac-se_frac,ymax=frac+se_frac)) + 
        ggplot2::geom_bar(stat="identity", color="black") + 
        ggplot2::geom_errorbar()
        )
}
