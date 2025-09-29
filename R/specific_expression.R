# Cleanup by Paul:
# moved check_input_se to utils
# removed alpha


#' computes specific expression and regulation
#' @description computes specific expression and regulation specific for certain groups of cells, weighting the individual cells accordingly.
#' normalization of the expression column is such that the total expression of a gene is approximately that of the bulk,
#' while the regulation column corresponds to the average expression of a single cell of that type.
#' This function provides the main computations accessed by specific_expression_regulation
#'
#' @param tissuemodel tissuemodel tibble as returned by fit_tissue
#' @param sclibrary single cell data matrix, columns are individual cells and must be indexable by the cell_ids contained in the celldf.
#' @param mapping dataframe that maps each cell (column name given by "cell_id") to a category (column name given by "by")
#' @param by name of the categorical variable defining groups of cells. defaults to "celltype"
#' @param weight name of the column that holds the fitted cell weights. defaults to "weight"
#' @param bulk_id name of the column that identifies the bulks. defaults to "bulk_id"
#' @param cell_id name of the column that contains the ids for the cells that match the column names of sclibrary. defaults to "cell_id"
#' @param compute_total if set to TRUE, also computes the total expression explained by ALL cell types. Can be used to compare resulting virtual tissue to the bulk.
#'
#' @return data frame holding the expression and regulation of each gene within each bulk and celltype.
#'
#' @examples
#' library(dplyr)
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
#' celldf <- fitted_tissue$tissuemodel %>% inner_join(mapping, by = "cell_id")
#' 
#' csre <- tissueResolver:::fast_specific_expression_regulation(
#'   tissuemodel = fitted_tissue,
#'   sclibrary = sc,
#'   mapping = mapping,
#'   by = "celltype",
#'   weight = "weight",
#'   bulk_id = "bulk_id",
#'   cell_id = "cell_id",
#'   compute_total = TRUE
#' )
fast_specific_expression_regulation <- function(

  tissuemodel,
  sclibrary,
  mapping,
  by = "celltype",
  weight = "weight",
  bulk_id = "bulk_id",
  cell_id = "cell_id",
  compute_total = TRUE 

) {
  # check if input is a tissuemodel
  if(!is_tissuemodel(tissuemodel)) {
    stop("Function called with an object that is not a tissuemodel")
  }

  # join the cell type column of mapping to tissuemodel tibble in order to comprise
  # the type of each cell
  celldf <- tissuemodel$tissuemodel %>% 
    dplyr::inner_join(mapping, by = cell_id)

  # transform sclibrary to matrix if it is a tibble
  if( tibble::is_tibble(sclibrary) ) {
    sclibrary <- df_to_matrix(sclibrary, cols = cell_id)
  }

  # check for correct formatting of resulting data frame celldf
  check_input_se(
    celldf, 
    sclibrary, 
    by, 
    weight, 
    bulk_id, 
    cell_id
    )

  # get genes of single cell library
  genes <- rownames(sclibrary)

  # for each cell type (by) and bulk (bulk_id) calculate the average expression 
  csre <- do.call(rbind,
    lapply(split(celldf, celldf[[by]]), function(df_ct) {

      do.call(rbind, 
        lapply(split(df_ct, df_ct[[bulk_id]]), function(df_ct_bulk) {
          
          # get current cell IDs and weights
          cellids <- df_ct_bulk[[cell_id]]
          weights <- stats::setNames(df_ct_bulk[[weight]], cellids)

          # specific expression = X_ctype beta_ctype
          se <- as.matrix(sclibrary[genes, cellids]) %*% as.matrix(weights[cellids])

          # bind specific expression and regulation for each bulk together
          tibble::tibble(
            bulk_id = df_ct_bulk[[bulk_id]] %>% unique(),
            by = df_ct_bulk[[by]] %>% unique(),
            gene = genes,
            expression = as.vector(se),
            regulation = as.vector(se / sum(weights)),
            weights = sum(weights)
            )

        })
      )
    })
    )

  # fill missing values in csre with 0
  csre <- csre %>%
    tidyr::complete(bulk_id, by, gene, 
      fill = list(expression = 0, regulation = 0, weights = 0)
      )

  # calculate the total (sum) expression and regulation for csre
  if(isTRUE(compute_total)) {
    csre_total <- csre %>% 
      dplyr::group_by(
        dplyr::across(c(bulk_id, "gene"))
        ) %>%
      dplyr::reframe(
        by = "total_explained",
        expression = sum(expression), 
        regulation = expression / sum(weights)
      )

    csre <- rbind(csre %>% dplyr::select(!weights), csre_total)
  }

  # format output to match input columns
  colnames(csre) <- c(
    bulk_id,
    by,
    colnames(csre)[3:length(colnames(csre))]
  )

  return(csre)
}

#' computes specific expression and regulation
#' @description computes specific expression and regulation specific for certain groups, weighting the individual cells accordingly.
#' normalization of the expression column is such that the total expression of a gene is approximately that of the bulk,
#' while the regulation column corresponds to the average expression of a single cell of that type.
#'
#' @param tissuemodel tissuemodel as returned by fit_tissue
#' @param sclibrary single cell data matrix, columns are individual cells and must be indexable by the cell_ids contained in the celldf.
#' @param mapping dataframe that maps each cell (column name given by cell_id) to a category (column name given by "by")
#' @param by name of the categorical variable defining groups of cells. defaults to "celltype"
#' @param weight name of the column that holds the fitted cell weights. defaults to "weight"
#' @param bulk_id name of the column that identifies the bulks. defaults to "bulk_id"
#' @param cell_id name of the column that contains the ids for the cells that match the column names of sclibrary. defaults to "cell_id"
#' @param compute_total if set to TRUE, also computes the total expression explained by ALL cell types. Can be used to compare resulting virtual tissue to the bulk.
#'
#' @return data frame holding the expression and regulation of each gene within each bulk and celltype.
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
#'   csre <- specific_expression_regulation(
#'     tissuemodel = fitted_tissue_boot,
#'     sclibrary = sc,
#'     mapping = mapping,
#'     by = "celltype",
#'     weight = "weight",
#'     bulk_id = "bulk_id",
#'     cell_id = "cell_id",
#'     compute_total = TRUE
#'   )
#'
#' @export
specific_expression_regulation <- function(

  tissuemodel,
  sclibrary,
  mapping,
  by = "celltype",
  weight = "weight",
  bulk_id = "bulk_id",
  cell_id = "cell_id",
  compute_total = TRUE 

) {

  # loop over each bootstrap run to calculate csre/wsce for each tm
  if("bootstrap" %in% names(tissuemodel) && tissuemodel$bootstrap == TRUE) {

    # get number of bootstrap runs
    n_boot <- tissuemodel$bootstrap_nruns

    # loop over all bootstrap runs and compute csre
    message(paste0("Computing CSRE for ", n_boot, " bootstrap runs: "))

    # set up progress bar
    pb <- utils::txtProgressBar(min = 0, max = n_boot, initial = 0, style = 3)

    for(run_i in 1:n_boot) {
      this_model <- list(
        "fitted_genes" = tissuemodel$fitted_genes, 
        "tissuemodel" = tissuemodel$tissuemodel[[run_i]]
        )

      # calculate the csre for the first bootstrap run
      curr_csre <- fast_specific_expression_regulation(
        this_model,
        sclibrary,
        mapping,
        by = by,
        weight = weight,
        bulk_id = bulk_id,
        cell_id = cell_id,
        compute_total = compute_total
      )

      if(run_i == 1) {

        csre <- curr_csre

      } else {

        # for all other iterations combine expression and regualtion as a list
        csre <- csre %>%
          dplyr::full_join(curr_csre, by = c(bulk_id, by, "gene")) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            expression = list(c(expression.x, expression.y)),
            regulation = list(c(regulation.x, regulation.y))
          ) %>%
          dplyr::select(-expression.x, -expression.y, -regulation.x, -regulation.y)

      }

      # update progress bar
      utils::setTxtProgressBar(pb, run_i)
    }

    # change NA values to 0 and add 0s for bootstrap runs where a cell 
    # type was not present that is present in later runs
    csre$expression <- lapply(csre$expression, function(row) {
        row[is.na(row)] <- 0
        c(rep(0, times = n_boot - length(row)), row)
      })

    csre$regulation <- lapply(csre$regulation, function(row) {
        row[is.na(row)] <- 0
        c(rep(0, times = n_boot - length(row)), row)
      })

    # calculate mean over all bootstrap runs
    csre <- csre %>%
      dplyr::rename(
        expression_boot = expression,
        regulation_boot = regulation
      ) %>%
      dplyr::mutate(
        expression = mean(expression_boot),
        regulation = mean(regulation_boot)
      ) %>%
      dplyr::ungroup()

      message("\t")
  } else {

    message("Computing CSRE for single tissue model")

    # if tissuemodel does not result from bootstrapping directly compute csre
    csre <- fast_specific_expression_regulation(
      tissuemodel,
      sclibrary,
      mapping,
      by = by,
      weight = weight,
      bulk_id = bulk_id,
      cell_id = cell_id,
      compute_total = compute_total
    )
  }

  return(csre)
}



#' computes TOTAL explained expression for quality control
#' @description computes the TOTAL explained expression, i.e., the specific gene expression within each bulk.
#' I.e., it uses the weights computed in fit tissue to determine the explained specific expression of each gene in each bulk sample without grouping cells.
#'
#' @param tissuemodel tissuemodel as returned by fit_tissue
#' @param sclibrary single cell data matrix, columns are individual cells and must be indexable by the cell_ids contained in the celldf.
#' @param by name of the categorical variable. Here it  is only used as index and the corresponding value in the returned data frame will always be "total explained".
#' In particular it is here not used to group cells. defaults to "celltype"
#' @param weight name of the column that holds the fitted cell weights. defaults to "weight"
#' @param bulk_id name of the column that identifies the bulks. defaults to "bulk_id"
#' @param cell_id name of the column that contains the ids for the cells that match the column names of sclibrary. defaults to "cell_id"
#'
#' @return data frame holding the explained expression for each gene in each bulk
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
#' csre_explained <- explained_expression(
#'             tissuemodel=fitted_tissue,
#'             sclibrary=sc,
#'             by="celltype",
#'             weight="weight",
#'             bulk_id = "bulk_id",
#'             cell_id="cell_id")
#'
#' @export
explained_expression <- function(
                          tissuemodel,
                          sclibrary,
                          by="celltype",
                          weight="weight",
                          bulk_id = "bulk_id",
                          cell_id="cell_id") {
  if(!is_tissuemodel(tissuemodel)) {
    stop("function called with an object that is not a tissuemodel")
  }
  if( tibble::is_tibble(sclibrary) ) {
    sclibrary <- df_to_matrix(sclibrary, cols=cell_id)
  }
  celldf <- tissuemodel$tissuemodel
  # do not check, if "by" is in celldf (it doesn't need to be)
  check_input_se(celldf, sclibrary, FALSE, weight, bulk_id, cell_id)

  # weigh each cell appearing in the tissuemodel with the corresponding weight computed by fit_tissue
  # then reshape to having the expression of each gene associated to each bulk as separate column
  # then group all entries which belong to the same gene AND the same bulk and compute the sum of expression (explained expression) over this group
  # finally cancel the weight and cell_id column, i.e., se contains bulk_id, gene and explained expression
  se <- cbind(
      celldf,
      celldf %>% dplyr::pull(weight) * t(sclibrary[,celldf %>% dplyr::pull(cell_id)])
      ) %>%
    tidyr::pivot_longer(
      -dplyr::all_of(c(bulk_id, cell_id, weight)),
      names_to="gene",
      values_to="expression"
      ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(bulk_id)), gene) %>%
    dplyr::mutate(expression = sum(expression)) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(c(bulk_id, "gene", "expression")))
  # add by column  with constant value total_explained       
  se[[by]] <- "total_explained"
  # discard all non unique entries
  return(
    se %>%
      dplyr::select(dplyr::all_of(c(bulk_id, by, "gene", "expression"))) %>%
      unique()
      )
}
