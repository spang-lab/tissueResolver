#' computes quality scores for each gene and for each bulk
#' @description Computes gene-specific and sample-specific quality scores:
#' This function computes the mean relative residual and avarage bootstrap variance for each gene and each bulk sample.
#'
#' @param csre the total explained expression obtained by fitting the sclibrary to the bulk tissue as given by explained_expression
#' @param bulks dataframe containing the gene expression of one bulk sample in each colum
#' 
#' @return a list of dataframes:
#' - genes: contains quality scores by gene + mean gene-wise actual and explained total expression
#' - bulks: contains quality scores to judge individual virtual tissues
#' 
#' @examples 
#' set.seed(1)
#' ngenes <- 4
#' nbulks <- 3
#' ncells <- 5
#' bulks <- matrix(runif(ngenes * nbulks), nrow = ngenes, ncol = nbulks)
#' sc <- matrix(runif(ngenes * ncells), nrow = ngenes, ncol = ncells)
#' bulknames <- paste0(rep("bulk_", nbulks), seq.int(1, nbulks))
#' cellnames <- paste0(rep("cell_", ncells), seq.int(1, ncells))
#' genenames_bulks <- paste0(rep("gene_", ngenes), seq.int(1, ngenes))
#' genenames_sc <- paste0(rep("gene_", ngenes), seq.int(1, ngenes))
#' rownames(bulks) <- genenames_bulks
#' rownames(sc) <- genenames_sc
#' colnames(bulks) <- bulknames
#' colnames(sc) <- cellnames
#' mapping <- tibble(
#'   cell_id = cellnames,
#'   celltype = c("CTA", "CTB", "CTA", "CTC", "CTA")
#' )
#' grouping <- tibble(
#'   "bulk_id" = c("bulk_1", "bulk_2", "bulk_3"),
#'   "group" = c("bgroup_1", "bgroup_2", "bgroup_1")
#' )
#' fitted_tissue <- fit_tissue(bulks, sc, bootstrap = TRUE, bootstrap_nruns = 5, bootstrap_pctcells = 50)
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
#' qc <- quality_scores(csre, bulks)
#'@export

quality_scores <- function(csre, bulks) {
  if( ! "total_explained" %in% (csre %>% pull(celltype) %>% unique())) {
    stop("need total explained gene expression in celldf. Call specific_expression_regulation with compute_total=TRUE.")
  }
  bulkdata <- promote_matrix(bulks)
  if( ncol(bulkdata) < 20 ) {
    warning("small number of bulks. quality scores will be unreliable.")
  }
  bulknames <- intersect(colnames(bulkdata), csre %>% pull(bulk_id) %>% unique())
  bulkdf <- as_tibble(bulkdata) %>% add_column(gene = rownames(bulkdata)) %>% pivot_longer(-gene, values_to = "bulk_expression", names_to = "bulk_id")
  total_fitted <- csre %>% filter(bulk_id %in% bulknames) %>% filter(celltype == "total_explained")  %>% ungroup() %>% dplyr::select(-celltype)
  total_fitted <- total_fitted %>% inner_join(bulkdf, by=c("gene", "bulk_id"))
  relres_score <- function(vec, bulkexpr) {
    # this function works both for vectors as well as scalars and can be used in the bootstrap and the non-bootstrap version.
    v = abs(vec - bulkexpr) / (vec + bulkexpr)
    v[is.na(v)] = 0.0 # if denom is zero, the fit is still perfect and there is no deviation.
    return(v)
  }
  if( "expression_boot" %in% names(csre) ) {
    # is a bootstrap run.
    res <- total_fitted %>% rowwise() %>%  
	    mutate(relres_boot = list(relres_score(expression_boot, bulk_expression))) %>% mutate(relres = mean(relres_boot)) %>% mutate(relres_var = var(relres_boot)) %>%  # relres
		ungroup() %>% dplyr::select(bulk_id, gene, relres, relres_var, expression, bulk_expression)
    gene_qc <- res %>%
      group_by(gene) %>%
      summarise(relres_sample_var = var(relres), relres = mean(relres), relres_mean_var = mean(relres_var), expression = mean(expression), bulk_expression = mean(bulk_expression))
    bulk_qc <- res %>% group_by(bulk_id) %>% summarise(relres_sample_var = var(relres), relres = mean(relres), relres_mean_var = mean(relres_var)) 
  } else {
    res <- total_fitted %>% mutate(relres = list(relres_score(expression, bulk_expression))) %>%
      ungroup() %>%
      dplyr::select(bulk_id, gene, relres, expression, bulk_expression) 

    gene_qc <- res %>% group_by(gene) %>% summarise(
                                            relres_sample_var = var(relres),
                                            relres = mean(relres),
                                            expression = mean(expression),
                                            bulk_expression=mean(bulk_expression))
    bulk_qc <- res %>% group_by(bulk_id) %>% summarise(
                                               relres_sample_var = var(relres),
                                               relres = mean(relres))
  }
  return(list(genes = gene_qc, bulks = bulk_qc))
}
