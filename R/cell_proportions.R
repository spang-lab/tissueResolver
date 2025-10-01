#' Compute cell type proportions from a single tissue model
#'
#' @description
#' Returns the sum over weights of certain cell populations (grouped by
#' mapping for each bulk) for a single tissue model retrieved by `fit_tissue()`
#' WITHOUT bootstrap.
#'
#' @param tissuemodel A `tissuemodel` as returned by `fit_tissue()` WITHOUT
#'    bootstrapping, or a data frame that contains cell ids, weights and a
#'    categorical variable, e.g., celltype.
#' @param mapping A dataframe that maps each cell (column name given by
#'    `cell_id`) to a category (column name given by `by`).
#' @param by The name of the categorical variable defining groups of cells.
#'    Defaults to `"celltype"`.
#' @param weight The name of the column that holds the fitted cell weights.
#'    Defaults to `"weight"`.
#' @param bulk_id The name of the column that identifies the bulks.
#'    Defaults to `"bulk_id"`.
#' @param cell_id The name of the column that identifies the cells.
#'    Defaults to `"cell_id"`.
#'
#' @return A data frame holding the cumulated weights for each cell population
#'    within each bulk.
#'
#' @examples
#' library(tibble)
#' set.seed(42)
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
cell_proportions_singlemodel <- function(
  tissuemodel,
  mapping = NULL,
  by = "celltype",
  weight = "weight",
  bulk_id = "bulk_id",
  cell_id = "cell_id"
) {

  if (is_tissuemodel(tissuemodel) && is.null(mapping)) {
    stop("tissuemodel must be a data frame if no mapping is given")
  } else if (is_tissuemodel(tissuemodel)) {
    check_input_mapping(tissuemodel, mapping, by, cell_id)

    # join the mapping from cell_id to cell type to the tissuemodel
    celldf <- tissuemodel$tissuemodel %>% 
      dplyr::inner_join(mapping, by = cell_id)
  } else { # is no tissuemodel (columns will be checked below)
    celldf <- tissuemodel
  }

  # input schecks
  for (col in c(by, weight, bulk_id, cell_id)) {
    if (!col %in% names(celldf)) {
      stop(paste0("Column (\"", col, "\") is not a name of celldf"))
    }
  }

  props <- celldf %>%
    # group cells declaring the same bulk and having the same by category
    # together
    dplyr::group_by(
      dplyr::across(dplyr::all_of(bulk_id)), dplyr::across(dplyr::all_of(by))
    ) %>%
    # sum the weights of the above groups
    dplyr::mutate(wsum = sum(weight)) %>%
    # discard redundant entries
    dplyr::select(dplyr::all_of(c(bulk_id, by, "wsum"))) %>%
    unique() %>%
    dplyr::rename(weight = wsum)

  return(props)

}

#' Compute cell type proportions from tissuemodel
#'
#' @description
#' Returns the sum over weights of certain cell populations (grouped by
#' mapping) for each bulk for a `tissuemodel` retrieved by `fit_tissue()`.
#'
#' @param tissuemodel A tissuemodel as returned by `fit_tissue()`, or a data
#'    frame that contains cell ids, weights and a categorical variable, e.g.,
#'    celltype.
#' @param mapping A dataframe that maps each cell (column name given by
#'    `cell_id`) to a category (column name given by `by`).
#' @param by The name of the categorical variable defining groups of cells.
#'    Defaults to `"celltype"`.
#' @param weight The name of the column that holds the fitted cell weights.
#'    Defaults to `"weight"`.
#' @param bulk_id The name of the column that identifies the bulks.
#'    Defaults to `"bulk_id"`.
#' @param cell_id The name of the column that identifies the cells.
#'    Defaults to `"cell_id"`.
#'
#' @return A data frame holding the cumulated weights for each cell population
#'    within each bulk for each bootstrap run as a list in the `weight_boot`
#'    column, if the tissuemodel included bootstrapping. The `weight` column
#'    holds the mean over the bootstrap runs across each bulk-celltype group
#'    and the last column holds the percentage of each cell type within a bulk.
#'
#' @examples
#' library(tibble)
#' set.seed(42)
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
#'
#' @export
#
cell_proportions <- function(
  tissuemodel,
  mapping = NULL,
  by = "celltype",
  weight = "weight",
  bulk_id = "bulk_id",
  cell_id = "cell_id"
) {

  if ("bootstrap" %in% names(tissuemodel) && tissuemodel$bootstrap == TRUE) {
    message(paste0(
      "Computing cell proportions for ",
      tissuemodel$bootstrap_nruns, " bootstrap runs:"
    ))

    # set up progress bar
    pb <- utils::txtProgressBar(
      min = 0,
      max = tissuemodel$bootstrap_nruns,
      initial = 0,
      style = 3
    )

    # get proportions for first bootstrap run
    irun <- 1

    # fetch the tissuemodel of the first bootstrap run
    this_model <- list(
      fitted_genes = tissuemodel$fitted_genes,
      tissuemodel = tissuemodel$tissuemodels[[irun]]
    )
    # compute the cumulated weights of this model for the cell population in
    # each bulk
    props <- cell_proportions_singlemodel(
      this_model,
      mapping,
      by,
      weight,
      bulk_id,
      cell_id
    ) %>%
      dplyr::rename(tmpweight_boot = weight)

    # loop over remaining bootstrap runs and concatenate values
    for (irun in 2:tissuemodel$bootstrap_nruns) {
      # fetch the tissuemodel of the irun-th bootstrap run
      this_model <- list(
        fitted_genes = tissuemodel$fitted_genes,
        tissuemodel = tissuemodel$tissuemodels[[irun]]
      )
      # compute the cumulated weights of this model for the cell population
      # in each bulk
      thisprops <- cell_proportions_singlemodel(
        this_model,
        mapping,
        by,
        weight,
        bulk_id,
        cell_id
      ) %>%
        dplyr::rename(current_weight = weight)
      # join the bootstrapped cell population weights as list to the props
      # dataframe within a column "tmpweight_boot"
      props <- props %>%
        dplyr::full_join(thisprops, by = c(bulk_id, by)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          tmpweight_boot = list(c(tmpweight_boot, current_weight)),
          # set NAs to 0, i.e. cells which are not drawn
          tmpweight_boot = list(tidyr::replace_na(tmpweight_boot, 0))
        ) %>%
        dplyr::select(-current_weight)

      # update progress bar
      utils::setTxtProgressBar(pb, irun)
    }

    # calcualte mean of weights across specific groups
    props <- props %>%
      dplyr::mutate(tmpweight = mean(tmpweight_boot)) %>%
      # dplyr::rename tmpweight with weight and tmp_boot with weight_boot
      dplyr::rename(
        !!weight := tmpweight,
        !!paste0(weight, "_boot") := tmpweight_boot
      )

    message("\t") # get to new line in terminal
  } else {
    message("Computing proportions for single tissue model")
    props <- cell_proportions_singlemodel(
      tissuemodel, mapping, by, weight, bulk_id, cell_id
    )
  }

  # calculate cell type proportions/fractions
  props <- props %>%
    dplyr::group_by(bulk_id) %>%
    dplyr::mutate(bulk_weight = sum(weight)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(per_cells = (weight / bulk_weight) * 100) %>%
    dplyr::select(!bulk_weight) %>%
    dplyr::rename(!!paste0(by, "_percent") := per_cells)

  return(props)

}
