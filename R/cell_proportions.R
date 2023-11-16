#' computes cell type proportions for single model (no bootstrap)
#' @description returns the sum over weights of certain cell populations (grouped by mapping for each bulk)
#' for a single tissue model retrieved by fit_tissue WITHOUT bootstrap.
#'
#' @param tissuemodel tissuemodel as returned by fit_tissue WITHOUT bootstrapping, or a data frame that contains cell ids, weights and a categorical variable, e.g., celltype
#' @param mapping dataframe that maps each cell (column name given by "cell_id") to a category (column name given by "by")
#' @param by name of the categorical variable defining groups of cells. defaults to "celltype"
#' @param weight name of the column that holds the fitted cell weights. defaults to "weight"
#' @param bulk_id name of the column that identifies the bulks. defaults to "bulk_id"
#' @param cell_id name of the column that identifies the cells. defaults to "cell_id"
#'
#' @return data frame holding the cumulated weights for each cell population within each bulk.
#' 
#' @examples 
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
#' fitted_tissue <- fit_tissue(bulks, sc)
#' 
#' mapping <- tibble(
#'   cell_id = cellnames,
#'   celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
#' )
#' 
#' cp <- tissueResolver:::cell_proportions_singlemodel(
#'   tissuemodel = fitted_tissue,
#'   mapping = mapping,
#'   by = "celltype",
#'   weight = "weight",
#'   bulk_id = "bulk_id",
#'   cell_id = "cell_id"
#' )
#'
#
cell_proportions_singlemodel <- function(
                                  tissuemodel,
                                  mapping = NULL,
                                  by="celltype",
                                  weight="weight",
                                  bulk_id = "bulk_id",
                                  cell_id="cell_id") {
  if( is_tissuemodel(tissuemodel) && is.null(mapping) ) {
    stop("tissuemodel must be a data frame if no mapping is given")
  } else if( is_tissuemodel(tissuemodel) ) {
    check_input_mapping(tissuemodel, mapping, by, cell_id)

    # join the mapping from cell_id to cell type to the tissuemodel
    celldf <- tissuemodel$tissuemodel %>% inner_join(mapping, by=cell_id)
  } else { # is no tissuemodel (columns will be checked below)
    celldf <- tissuemodel
  }
  if( ! by %in% names(celldf) ) {
    stop(paste0("by-selection (\"", by, "\") is not a name of celldf"))
  }
  if( ! weight %in% names(celldf) ) {
    stop(paste0("weights (\"", weights, "\") is not a name of celldf"))
  }
  if( ! bulk_id %in% names(celldf) ) {
    stop(paste0("bulk_id (\"", bulk_id, "\") is not a name of celldf"))
  }

  props <- celldf %>%
    # group cells declaring the same bulk and having the same by category together
    dplyr::group_by(dplyr::across(dplyr::all_of(bulk_id)), dplyr::across(dplyr::all_of(by))) %>%
    # sum the weights of the above groups
    mutate(wsum=sum(weight)) %>% 
    # discard redundant entries
    dplyr::select(dplyr::all_of(c(bulk_id, by, "wsum"))) %>% unique() %>%
    rename(weight = wsum)
  return(props)
}

#' computes cell type proportions from virtual tissues and label mappings
#' @description returns the sum over weights of certain cell populations for each bootstrap (if enabled) run performed in fit_tissue
#'
#' @param tissuemodel tissuemodel as returned by fit_tissue, or a data frame that contains cell ids, weights and a categorical variable, e.g., celltype
#' @param mapping dataframe that maps each cell (column name given by "cell_id") to a category (column name given by "by")
#' @param by name of the categorical variable defining groups of cells. defaults to "celltype"
#' @param weight name of the column that holds the fitted cell weights. defaults to "weight"
#' @param bulk_id name of the column that identifies the bulks. defaults to "bulk_id"
#' @param cell_id name of the column that identifies the cells. defaults to "cell_id"
#'
#' @return data frame holding the cumulated weights for each cell population within each bulk for each bootstrap run as a list in the weight_boot column.
#' The weight column holds the mean over the bootstrap runs across each bulk-celltype group.
#'
#' @examples
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
#' fitted_tissue_boot <- fit_tissue(
#'   bulks,
#'   sc,
#'   bootstrap = TRUE,
#'   bootstrap_nruns = 5,
#'   bootstrap_pctcells = 50
#' )
#' 
#' 
#' 
#' mapping <- tibble(
#'   cell_id = cellnames,
#'   celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
#' )
#' 
#' cell_proportions_boot <- cell_proportions(
#'   tissuemodel = fitted_tissue_boot,
#'   mapping = mapping,
#'   by = "celltype",
#'   weight = "weight",
#'   bulk_id = "bulk_id",
#'   cell_id = "cell_id"
#' )
#' @export
#
cell_proportions <- function(
                      tissuemodel,
                      mapping = NULL,
                      by="celltype",
                      weight="weight",
                      bulk_id = "bulk_id",
                      cell_id="cell_id") {
  if ("bootstrap" %in% names(tissuemodel) && tissuemodel$bootstrap == TRUE) {
    irun <- 1
    message(paste0("Computing cell proportions for bootstrap sample ", irun, " / ", tissuemodel$bootstrap_nruns))
    # fetch the tissuemodel of the first bootstrap run
    this_model <- list(
                    fitted_genes = tissuemodel$fitted_genes,
                    tissuemodel = tissuemodel$tissuemodels[[irun]])
    # compute the cumulated weights of this model for the cell population in each bulk
    props <- cell_proportions_singlemodel(
                            this_model,
                            mapping,
                            by,
                            weight,
                            bulk_id, 
                            cell_id) %>% rename (tmpweight_boot = weight)
    for( irun in 2:tissuemodel$bootstrap_nruns ){
      message(paste0("Computing cell proportions for bootstrap sample ", irun, " / ", tissuemodel$bootstrap_nruns))
      # fetch the tissuemodel of the irun-th bootstrap run
      this_model <- list(
                      fitted_genes = tissuemodel$fitted_genes,
                      tissuemodel = tissuemodel$tissuemodels[[irun]])
      # compute the cumulated weights of this model for the cell population in each bulk
      thisprops <- cell_proportions_singlemodel(
                      this_model,
                      mapping,
                      by,
                      weight,
                      bulk_id,
                      cell_id) %>% rename(current_weight = weight)
      # join the bootstrapped cell population weights as list to the props dataframe
      # within a column "tmpweight_boot"
      props <-  props %>%
                full_join(thisprops, by=c(bulk_id, by)) %>%
                rowwise() %>%
                mutate(tmpweight_boot = list(c(tmpweight_boot, current_weight))) %>%
                # set NAs to 0, i.e. cells which are not drawn
                mutate(tmpweight_boot = list(replace_na(tmpweight_boot,0))) %>%
                dplyr::select(-current_weight)

    }

    props <-  props %>% 
              mutate(tmpweight = mean(tmpweight_boot)) %>%
              # rename tmpweight with weight and tmp_boot with weight_boot
              rename(!!weight := tmpweight, !!paste0(weight, "_boot") := tmpweight_boot)
              # quoname is deprecated
              # rename(!!quo_name(weight) := tmpweight) %>%
              # rename(!!quo_name(paste0(weight, "_boot")) := tmpweight_boot)
    return(props)
  } else {
    return(cell_proportions_singlemodel(tissuemodel, mapping, by, weight, bulk_id, cell_id))
  }
}
