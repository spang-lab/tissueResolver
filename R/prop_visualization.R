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
    if (!groupA %in% (grouping %>% pull(group) %>% unique()) ||
        !groupB %in% (grouping %>% pull(group) %>% unique())) {
        stop("given groups are not levels in the grouping file.")
    }
    if( ! "group" %in% names(grouping) ) {
      stop("no group in grouping df")
    }
    if( ! "bulk_id" %in% names(grouping) ) {
      stop("no bulk_id in grouping df")
    }


    groupvec <- grouping %>% pull(group)
    names(groupvec) <- grouping %>% pull(bulk_id)
    if( "weight_boot" %in% names(props) ){
        propdf <- props %>% drop_na() %>%
            dplyr::select(-weight_boot)
    } else {
        propdf <- props %>% drop_na()
    }
    if( ! "group" %in% names(props)) {
        propdf <- propdf %>% inner_join(grouping, by="bulk_id")
    }
    propdf <- propdf %>%
        pivot_wider(names_from = bulk_id, values_from=weight, values_fill = 0.0) %>%
        pivot_longer(c(-celltype, -group), values_to = "weight", names_to="bulk_id") %>%
        filter(groupvec[bulk_id] == group) %>% group_by(group, celltype) %>%
        summarise(mweight = mean(weight), se=sd(weight)/sqrt(length(weight))) %>%
        pivot_wider(names_from = group, values_from=c(mweight, se)) %>%
        mutate(frac = .data[[paste0("mweight_", groupA)]] / .data[[paste0("mweight_", groupB)]],
               se_frac = .data[[paste0("se_", groupA)]] / .data[[paste0("mweight_", groupA)]] +
                 .data[[paste0("se_", groupB)]] / .data[[paste0("mweight_", groupB)]])
    return(ggplot(propdf, aes(x=celltype, y=frac, ymin=frac-se_frac,ymax=frac+se_frac)) + geom_bar(stat="identity", color="black") + geom_errorbar())
}
