#' plots cell type specific expression comparing two bulk groups
#'
#' @description This function generates several bar and heatmapplots for comparing cell type specific expression comparing two bulk groups.
#' 
#' @param csre object returned by specific_expression_regulation, with an added bulk group column
#' @param heatmapviz Quantity to show in the heat map.
#'             May be either "logdiffexpr" (for the logarithmic differential expression),
#'             "logdiffreg" (log differential regulation, per cell), or
#'             "relchange" (relative change in expression compared to the average of the two groups), or
#'             "relregchange" (relative change in regulation -- per cell -- compared to the average of the two groups)
#' @param barplotviz transformation to show in the bar plot for each gene
#'             May be either 
#'             "linear", which will show the reconstructed gene expression of the two bulk groups as a stacked colored bar plot, or
#'             "log", showing the total log-expression in both groups, or
#'             "logdiff", showing the differential log-expression between the two groups in a single bar plot,
#'             "relative", or "relativeA" showing the group-A expression as 100% fractions and the group-B relative to that.
#'             "relativeB", similarly, but normalized to the group-B expression
#'
#' @param groupA identifier of bulk group A (factor in the group column)
#' @param groupB identifier of bulk group B (factor in the group column)
#' @param ctypes if non-NULL, the order in which the cell types will be displayed.
#'             All celltypes present in the dataset but not in this array will be summed as "other".
#' @param addexpressionlevel if set to TRUE and barplotviz == relative, add an additional bar plot that indicates the expression level of the gene present in bulk group A
#'
#' @return ggplot object
#' 
#' @examples 
#' library(gridExtra)
#' set.seed(1)
#' ngenes <- 4
#' nbulks <- 3
#' ncells <- 5
#' 
#' bulks <- matrix(runif(ngenes * nbulks), nrow = ngenes, ncol = nbulks)
#' sc <- matrix(runif(ngenes * ncells), nrow = ngenes, ncol = ncells)
#' bulknames <- paste0(rep("bulk_", nbulks), seq.int(1, nbulks))
#' cellnames <- paste0(rep("cell_", ncells), seq.int(1, ncells))
#' genenames_bulks <- paste0(rep("gene_", ngenes), seq.int(1, ngenes))
#' genenames_sc <- paste0(rep("gene_", ngenes), seq.int(1, ngenes))
#' 
#' rownames(bulks) <- genenames_bulks
#' rownames(sc) <- genenames_sc
#' colnames(bulks) <- bulknames
#' colnames(sc) <- cellnames
#' 
#' mapping <- tibble(
#'   cell_id = cellnames,
#'   celltype = c("CTA", "CTB", "CTA", "CTC", "CTA")
#' )
#' 
#' grouping <- tibble(
#'   "bulk_id" = c("bulk_1", "bulk_2", "bulk_3"),
#'   "group" = c("bgroup_1", "bgroup_2", "bgroup_1")
#' )
#' 
#' fitted_tissue <- fit_tissue(
#'   bulks,
#'   sc
#' )
#' 
#' csre <- specific_expression_regulation(
#'   tissuemodel = fitted_tissue,
#'   sclibrary = sc,
#'   mapping = mapping,
#'   by = "celltype",
#'   weight = "weight",
#'   bulk_id = "bulk_id",
#'   cell_id = "cell_id",
#'   compute_total = TRUE
#' )
#' 
#' csre <- csre %>%
#'   inner_join(grouping %>% dplyr::select(bulk_id, group)) %>%
#'   drop_na()
#' 
#' 
#' csre_plot <- plot_csre(
#'   csre = csre,
#'   heatmapviz = "relchange",
#'   barplotviz = "relative",
#'   groupA = "bgroup_1",
#'   groupB = "bgroup_2",
#'   ctypes = c("CTA", "CTC"),
#'   addexpressionlevel = TRUE
#' )
#' 
#' print(csre_plot)
#' @export
#'
plot_csre <- function(
                csre,
                heatmapviz = "relchange",
                barplotviz="linear",
                groupA = NA,
                groupB = NA,
                ctypes = NULL,
                addexpressionlevel = FALSE ) {
  
  # retrieve bulk groups from csre tibble (e.g. GCB vs. ABC)
  groups <- csre %>% pull(group) %>% unique() %>% sort()
  if( is.na(groupA) || is.na(groupB) ){
    stop("no groups given.")
  }

  if( ! is.null(ctypes) ) {
    # remove_celltypes will store the "other" category of cells, i.e.,
    # cells which do not belong to the ctypes array
    remove_celltypes <- setdiff(
      csre %>%
        filter(celltype != "total_explained") %>%
        pull(celltype) %>% unique(),
      ctypes)
  } else {
    remove_celltypes <- c()
  }

  # average over all bulks within the same (group, celltype, gene)-group:
  groupavg <- csre %>%
    # assign to the celltypes belonging to remove_celltypes the value "other"
    mutate(celltype = ifelse(celltype %in% remove_celltypes, "other", celltype)) %>%
    # set regulation entries which are NaN to 0
    mutate(regulation = ifelse(is.na(regulation), 0.0, regulation)) %>%
    # set expression entries which are NaN to 0
    mutate(expression = ifelse(is.na(expression), 0.0, expression)) %>%
    # group across bulk_id, celltype, gene and group
    group_by(bulk_id, celltype, gene, group) %>%
    # compute cumulated sums over these groups
    mutate(expression = sum(expression)) %>%
    mutate(regulation = sum(regulation)) %>%
    unique() %>% ungroup() %>%
    # consider only the two bulk groups A and B
    filter(group %in% c(groupA, groupB)) %>%
    # group across bulk group, gene and celltype
    group_by(group,gene,celltype) %>%
    # compute means across this group,
    # i.e., compare bulks belonging to same bulkgroup
    mutate(avgexpr = mean(expression)) %>%
    mutate(avgregulation = mean(regulation)) %>%
    ungroup() %>% dplyr::select(group, gene, celltype, avgexpr, avgregulation) %>% unique() %>% drop_na()

  # the average total explained expression across bulks for each (group,gene)-group
  # this is just the lower part of groupavg but without the regulation column
  total_explained <- groupavg %>% filter(celltype == "total_explained") %>% dplyr::select(group, gene, avgexpr)


  diff <- groupavg %>%
    # join group avarage and total explained avarage 
    # and add the total_expr column giving the total expression across cells
    inner_join(total_explained %>% rename(total_expr = avgexpr), by=c("group", "gene")) %>%
    # give avgexpr, avgreg and total_expr separate columns for both bulkgroups A and B
    pivot_wider(names_from = group, values_from=c(avgexpr, avgregulation, total_expr), names_prefix="group_", values_fill = 0.0) %>%
    rename(avgexpr_A = paste0("avgexpr_group_", groupA)) %>%
    rename(avgexpr_B = paste0("avgexpr_group_", groupB)) %>%
    rename(avgregulation_A = paste0("avgregulation_group_", groupA)) %>%
    rename(avgregulation_B = paste0("avgregulation_group_", groupB)) %>%
    rename(total_expr_A = paste0("total_expr_group_", groupA)) %>%
    rename(total_expr_B = paste0("total_expr_group_", groupB)) %>%
    # for each row (comprising gene,celltype) give logp1 difference of avgexpr, avgreg and total_expr
    mutate(logdiffexpr = log2((avgexpr_A+1.0) / (avgexpr_B+1.0))) %>%
    mutate(logdiffreg = log2((avgregulation_A+1.0) / (avgregulation_B+1.0))) %>%
    mutate(logdifftotal = log2((total_expr_A+1.0)/(total_expr_B+1.0))) %>%
    # give relative difference between avgexpr relative to the tota_expr
    mutate(relchange = 0.5*(avgexpr_A-avgexpr_B)/(total_expr_A+total_expr_B)) %>%
    # give relative difference between avgregulation
    mutate(relregchange = 0.5*(avgregulation_A-avgregulation_B)/(avgregulation_A+avgregulation_B)) %>%
    # replace NAs by 0
    mutate(relregchange = ifelse(is.na(relregchange), 0.0, relregchange))

  # for the column specified in heatmapviz (e.g., "relchange") we compute the absolute value 5% and 95% quantile
  # then we take the maximum of both values and scale it with 1.1
  # this quantity will be used for heatmap plotting at the end
  limits <- 1.1*max(abs(diff %>% pull(heatmapviz) %>% quantile(probs=c(0.05, 0.95), na.rm = TRUE)))

  # we cluster based on genes and celltypes:
  tmp <- diff %>% dplyr::select(celltype, gene, all_of(heatmapviz)) %>% drop_na() %>%
    # give each gene a separate column storing the heatmapviz value for each celltype
    pivot_wider(celltype, names_from=gene, values_from=all_of(heatmapviz)) %>%
    # only consider the actual celltypes so neither "other" nor "total_explained"
    filter(! celltype %in% c("other", "total_explained"))
  # format tmp as matrix with celltypes as rownames  
  m <- tmp %>% dplyr::select(-celltype) %>% as.matrix()
  rownames(m) <- tmp %>% pull(celltype)

  # perform hierarchical clustering in order retrieve
  # hierarchical order of genes
  ord_genes  <- colnames(m)[hclust(dist(t(m)), method="median")$order]

  if( is.null(ctypes) ) {
    # hierarchically order ctypes
    ctypes <- rownames(m)[hclust(dist(m), method="median")$order]
    if( "other" %in% (diff %>% pull(celltype) %>% unique()) ) {
      # add "other" to the end:
      ctypes <- c(ctypes, "other")
    }
  } else {
    # if ctypes is not null keep the ordering of ctypes and add "other" to the end
    if( ! ("other" %in% ctypes) && ("other" %in% (diff %>% pull(celltype) %>% unique())) ) {
      ctypes <- c(ctypes, "other")
    }
  }


  # preparing the expression bar plots,
  # according to the user specified visualization
  if( barplotviz == "linear" ) {
    groupavg <- groupavg %>% filter(celltype != "total_explained") %>%
      # encode group column as factor type
      mutate(across(group, factor, levels=c(groupA, groupB)))
    exprchanges1 <-
      groupavg %>%
      # create for each bulk group a barplot, where for each gene the avgexpression contribution of each cell type
      # is visualized by the filling of the bar
      ggplot(aes(x=avgexpr, y=factor(gene, rev(ord_genes)), fill=factor(celltype, levels = ctypes))) +
      geom_bar(stat="identity", position="stack") +
      facet_wrap(~group, strip.position = "left") +
      ylab(NULL) + xlab(NULL) +
      theme(legend.position = "none")
  } else if( barplotviz == "log" ) {
    # logp1 transform the avg total expression
    plotdf <- total_explained %>%
      mutate(logexpr = log2(avgexpr+1.0)) %>% dplyr::select(-avgexpr) %>%
      mutate(across(group, factor, levels=c(groupA, groupB)))
    # create for each bulk group a barplot, where for each gene the logp1 avg total expression
    # is visualized by the length of the bar
    exprchanges1 <- plotdf %>%
      ggplot(aes(x=logexpr, y=factor(gene, rev(ord_genes)))) +
      geom_bar(stat="identity", position="stack") +
      facet_wrap(~group, strip.position = "left") +
      ylab(NULL) + xlab(NULL) +
      theme(legend.position = "none")
  } else if( barplotviz == "logdiff") {
    # for the total avarage expression plot
    # the log1p difference between groups A and B for each gene
    # the difference is visualized by the length of the bar
    exprchanges1 <-
      total_explained %>%
      mutate(logexpr = log2(avgexpr+1.0)) %>%
      dplyr::select(-avgexpr) %>%
      pivot_wider(names_from=group, values_from=logexpr) %>%
      rename(group_A = groupA, group_B = groupB) %>%
      mutate(logdiff=group_A-group_B) %>%
      ggplot(aes(x=logdiff, y=factor(gene, rev(ord_genes)))) +
      geom_bar(stat="identity", position="stack") +
      ylab(NULL) + xlab(NULL) +
      theme(legend.position = "none")

  } else if( barplotviz == "relativeA" || barplotviz == "relativeB" || barplotviz == "relative") {
    if( barplotviz == "relativeB" ) {
      relgroup = "group_B"
    } else {
      relgroup = "group_A"
    }
    # this dataframe stores for each gene within each bulkgroup the
    # fraction between the avgexpr of each celltype relative to the total expression
    # cumulated across celltypes within the bulkgroup specified in barplotviz ("relative" defaults to "relativeA")
    plotdf <-
      groupavg %>% dplyr::select(-avgregulation) %>%
      filter(celltype != "total_explained") %>%
      # give both groups their own column
      pivot_wider(names_from = group, values_from =  avgexpr, values_fill=0) %>%
      rename(group_A = all_of(groupA), group_B = all_of(groupB)) %>%
      group_by(gene) %>%
        # cumulate the avgexpr for each gene in the relgroup across celltypes
        mutate(normsum = sum(.data[[relgroup]])) %>%
        # frac gives the proportion between the avgexpr for each celltype
        # relative to the whole expression
        mutate(frac_A = group_A / normsum) %>%
        mutate(frac_B = group_B / normsum) %>%
      ungroup() %>% dplyr::select(-group_A, -group_B, -normsum) %>%
      # store only the frac values
      pivot_longer(c(frac_A, frac_B), names_prefix="frac_", names_to="group", values_to="frac")
    # give both groups their original name
    plotdf[,"group"] <- sapply(
      plotdf %>% pull(group),
      function(x){if(x == "A"){return(groupA)}else if(x=="B"){return(groupB)}else{return(NA)}}
    )

    # encode group column as factor type
    plotdf <- plotdf %>% mutate(group = factor(group, levels=c(groupA, groupB)))
    # as specified in the barplotpviz parameter the fractions are relative to one of both groups, justifying the
    # equal length of one of the barplots
    exprchanges1 <- plotdf %>%
      ggplot(aes(x=frac, y=factor(gene, rev(ord_genes)), fill=factor(celltype, levels = rev(ctypes)))) +
      geom_bar(stat="identity", position="stack") +
      facet_wrap(~group, strip.position = "left") +
      ylab(NULL) + xlab(NULL) +
      theme(legend.position = "none")

    if( addexpressionlevel ) {
      # consider only the avgexpr of groupA
      # logp1 transform it
      exprlevelplot <-
        total_explained %>%
        filter(group == all_of(groupA)) %>%
        dplyr::select(-group) %>%
        mutate(logexpr = log2(avgexpr+1)) %>%
        # the barplot show the logp1 transformed avgexpression of groupA
        ggplot(aes(x=logexpr, y=factor(gene, rev(ord_genes)))) +
        geom_bar(stat="identity", position="stack") +
        ylab(NULL) + xlab(NULL) +
        theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

      #join both plots
      exprchanges1 <- grid.arrange(exprchanges1,
                                   exprlevelplot,
                                   ncol=2, widths=c(0.8, 0.2))
    }
  } else {
    error("unknown key in barplotviz")
  }
  # plot the relchange of each gene for each celltype as tile plot
  # each tile corresponds to one gene and the color expresses the relchange stored in diff
  # relchange is the relative change between both bulkgroups in avgexpr
  exprchanges2 <- diff %>%
    filter(celltype != "total_explained") %>%
    ggplot(aes(y=factor(gene, rev(ord_genes)), x = factor(celltype, levels = ctypes), fill=.data[[heatmapviz]])) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-limits, limits), oob=scales::squish) +
    ylab(NULL) + xlab(NULL) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", plot.margin=margin(0,0,0,2,"cm"))

  p <- grid.arrange(exprchanges1, exprchanges2, ncol=2, widths=c(0.5, 0.5))

  return(p)
}
