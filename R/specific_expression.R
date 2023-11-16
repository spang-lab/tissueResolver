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
                                          by="celltype",
                                          weight="weight",
                                          bulk_id = "bulk_id",
                                          cell_id="cell_id",
                                          compute_total = TRUE) {
  if(!is_tissuemodel(tissuemodel)) {
    stop("function called with an object that is not a tissuemodel")
  }
  # check for correct formatting of mapping
  check_input_mapping(tissuemodel, mapping, by, cell_id)
  if( is_tibble(sclibrary) ) {
    sclibrary <- df_to_matrix(sclibrary, cols=cell_id)
  }
  # join the by column of mapping to tissuemodel tibble in order to comprise
  # the type of each cell
  celldf <- tissuemodel$tissuemodel %>% inner_join(mapping, by=cell_id)

  # check for correct formatting of resulting data frame celldf
  check_input_se(celldf, sclibrary, by, weight, bulk_id, cell_id)

  genes <- rownames(sclibrary)

  res <- tibble()
  # store all possible celltypes determined by "by" 
  ctypes <- celldf %>% pull(by) %>% unique()
  if( compute_total ) {
    # allow for computation of total expression explained by ALL celltypes
    ctypes <- c(ctypes, "total_explained")
  }
  # for each celltype loop through each bulk to retrive
  # subframe of celldf (called thesecells) holding only the specified celltype(s)
  for(ctype in ctypes) {
    for(bulkid in celldf %>% pull(bulk_id) %>% unique()) {
      # filter ALL cells with matching bulk_id
      if( ! is.na(ctype) && ctype == "total_explained" ){
        thesecells <- celldf %>% filter(across(bulk_id) == bulkid)
        # filter ONLY cells with matching celltype and matching bulk_id
      } else if( ! is.na(ctype) ) {
        thesecells <- celldf %>% filter(across(by) == ctype) %>% filter(across(bulk_id) == bulkid)
        # treat NA as celltype to filter by
      } else {
        thesecells <- celldf %>% filter(is.na(across(by)) ) %>% filter(across(bulk_id) == bulkid)
      }
      weights <- thesecells %>% pull(weight)
      cellids <- thesecells %>% pull(cell_id)
      names(weights) <- cellids
      # specific expression = X_ctype beta_ctype
      se <- as.matrix(sclibrary[genes,cellids]) %*% as.matrix(weights[cellids])
      # bind specific expression and regulation for each bulk together
      res <- res %>% rbind(tibble(
                       bulk_id = bulkid,
                       ctype = ctype,
                       gene = genes,
                       expression  = as.vector(se),
                       regulation = as.vector(se / sum(weights))
                     ))
    }
  }
  # just some formatting for clearer structure
  namesres <- names(res)
  namesres <- c(bulk_id, by, namesres[3:length(namesres)])
  names(res) <- namesres
  return(res)
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
                                    by="celltype",
                                    weight="weight",
                                    bulk_id = "bulk_id",
                                    cell_id="cell_id",
                                    compute_total = TRUE
                                  ) {
  if("bootstrap" %in% names(tissuemodel) && tissuemodel$bootstrap == TRUE) {

    message(paste0("Computing CSRE for bootstrap sample ", 1, " / ", tissuemodel$bootstrap_nruns))

    # retrieve model for each bootstrap run and compute csre
    this_model <- list(
                    fitted_genes = tissuemodel$fitted_genes,
                    tissuemodel = tissuemodel$tissuemodels[[1]])
    csre <- fast_specific_expression_regulation(
                    this_model,
                    sclibrary,
                    mapping,
                    by=by,
                    weight=weight,
                    bulk_id = bulk_id,
                    cell_id = cell_id,
                    compute_total = compute_total)

    for( irun in 2:tissuemodel$bootstrap_nruns ){
      message(paste0("Computing CSRE for bootstrap sample ", irun, " / ", tissuemodel$bootstrap_nruns))

      this_model <- list(
                      fitted_genes = tissuemodel$fitted_genes,
                      tissuemodel = tissuemodel$tissuemodels[[irun]])
      this_csre <- fast_specific_expression_regulation(
                      this_model,
                      sclibrary,
                      mapping,
                      by=by,
                      weight=weight,
                      bulk_id = bulk_id,
                      cell_id = cell_id,
                      compute_total = compute_total)
      # join csre's together and store resulting expression and regulation as list
      csre <- csre %>% right_join(this_csre, by =c(bulk_id, by, "gene")) %>%
        rowwise(all_of(c(bulk_id, by, "gene"))) %>%
        mutate(expression = list(c(expression.x, expression.y))) %>%
        mutate(regulation = list(c(regulation.x, regulation.y))) %>%
        dplyr::select(-expression.x, -expression.y, -regulation.x, -regulation.y)
    }
    # add column storing mean over all bootstrapped expressions and regulations
    csre <- csre %>% rename(expression_boot = expression) %>%
      mutate(expression = mean(expression_boot)) %>%
      rename(regulation_boot = regulation) %>%
      mutate(regulation = mean(regulation_boot))
    return(csre)
  } else {
    # if tissuemodel does not result from bootstrapping directly compute csre
    return(fast_specific_expression_regulation(tissuemodel, sclibrary, mapping, by=by, weight=weight, bulk_id = bulk_id, cell_id = cell_id, compute_total = compute_total))
  }
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
  if( is_tibble(sclibrary) ) {
    sclibrary <- df_to_matrix(sclibrary, cols=cell_id)
  }
  celldf <- tissuemodel$tissuemodel
  # do not check, if "by" is in celldf (it doesn't need to be)
  check_input_se(celldf, sclibrary, FALSE, weight, bulk_id, cell_id)

  # weigh each cell appearing in the tissuemodel with the corresponding weight computed by fit_tissue
  # then reshape to having the expression of each gene associated to each bulk as separate column
  # then group all entries which belong to the same gene AND the same bulk and compute the sum of expression (explained expression) over this group
  # finally cancel the weight and cell_id column, i.e., se contains bulk_id, gene and explained expression
  se <- cbind(celldf,
              celldf %>% pull(weight) * t(sclibrary[,celldf %>% pull(cell_id)])
              ) %>%
            pivot_longer(-dplyr::all_of(c(bulk_id, cell_id, weight)), names_to="gene", values_to="expression") %>%
            group_by(across(dplyr::all_of(bulk_id)), gene) %>% 
            mutate(expression = sum(expression)) %>%
            ungroup() %>%
            dplyr::select(dplyr::all_of(c(bulk_id, "gene", "expression")))
  # add by column  with constant value total_explained       
  se[[by]] <- "total_explained"
  # discard all non unique entries
  return(se %>% dplyr::select(dplyr::all_of(c(bulk_id, by, "gene", "expression"))) %>% unique())
}
