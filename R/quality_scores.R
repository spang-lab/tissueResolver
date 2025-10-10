#' Compute quality scores for genes and bulks
#'
#' @description
#' This function computes gene-specific and sample-specific quality scores.
#' Specifically, it computes the mean relative residual and average bootstrap
#' variance for each gene and each bulk sample.
#'
#' @param csre The cell type specific expression (`csre`) obtained by
#'    `specific_expression_regulation()`. Note that `csre` must include the
#'    total explained gene expression to calculate the quality scores.
#' @param bulks A matrix or dataframe containing the gene expression of bulk
#'    samples in each column and genes in each row. Row names must match the
#'   gene names in `csre`.
#'
#' @return A list of dataframes:
#' - `genes`: Contains quality scores by gene + mean gene-wise actual and
#'    explained total expression
#' - `bulks`: Contains quality scores to judge individual virtual tissues
#'
#' @examples
#' library(tibble)
#' set.seed(42)
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
#' fitted_tissue <- fit_tissue(
#'   bulks,
#'   sc,
#'   bootstrap = TRUE,
#'   bootstrap_nruns = 5,
#'   bootstrap_pctcells = 50
#' )
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
#'
#' @export
#'
quality_scores <- function(
  csre,
  bulks
) {

  if (!"total_explained" %in% (csre %>% dplyr::pull(celltype) %>% unique())) {
    stop(
      "need total explained gene expression in celldf. 
      Call specific_expression_regulation with compute_total=TRUE."
    )
  }

  bulkdata <- promote_matrix(bulks)

  if (ncol(bulkdata) < 20) {
    warning("small number of bulks. quality scores will be unreliable.")
  }

  bulknames <- intersect(
    colnames(bulkdata),
    csre %>% dplyr::pull(bulk_id) %>% unique()
  )

  bulkdf <- tidyr::as_tibble(bulkdata) %>%
    tibble::add_column(gene = rownames(bulkdata)) %>%
    tidyr::pivot_longer(
      -gene,
      values_to = "bulk_expression",
      names_to = "bulk_id"
    )

  total_fitted <- csre %>%
    dplyr::filter(bulk_id %in% bulknames) %>%
    dplyr::filter(celltype == "total_explained") %>%
    dplyr::ungroup() %>%
    dplyr::select(-celltype)

  total_fitted <- total_fitted %>%
    dplyr::inner_join(bulkdf, by = c("gene", "bulk_id"))

  relres_score <- function(vec, bulkexpr) {
    # this function works both for vectors as well as
    # scalars and can be used in the bootstrap and the non-bootstrap version.
    v <- abs(vec - bulkexpr) / (vec + bulkexpr)
    # if denom is zero, the fit is still perfect w/o deviation
    v[is.na(v)] <- 0.0
    v
  }

  if ("expression_boot" %in% names(csre)) {

    # is a bootstrap run.
    res <- total_fitted %>%
      dplyr::rowwise() %>%
      dplyr::mutate(relres_boot = list(
        relres_score(expression_boot, bulk_expression)
      )) %>%
      dplyr::mutate(relres = mean(relres_boot)) %>%
      dplyr::mutate(relres_var = stats::var(relres_boot)) %>%  # relres
      dplyr::ungroup() %>%
      dplyr::select(
        bulk_id, gene, relres, relres_var, expression, bulk_expression
      )

    gene_qc <- res %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        relres_sample_var = stats::var(relres),
        relres = mean(relres),
        relres_mean_var = mean(relres_var),
        expression = mean(expression),
        bulk_expression = mean(bulk_expression)
      )

    bulk_qc <- res %>%
      dplyr::group_by(bulk_id) %>%
      dplyr::summarise(
        relres_sample_var = stats::var(relres),
        relres = mean(relres),
        relres_mean_var = mean(relres_var)
      )

  } else {

    # no bootstrap run
    res <- total_fitted %>%
      dplyr::mutate(
        relres = relres_score(expression, bulk_expression)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        bulk_id, gene, relres, expression, bulk_expression
      )

    gene_qc <- res %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        relres_sample_var = stats::var(relres),
        relres = mean(relres),
        expression = mean(expression),
        bulk_expression = mean(bulk_expression)
      )

    bulk_qc <- res %>% 
      dplyr::group_by(bulk_id) %>%
      dplyr::summarise(
        relres_sample_var = stats::var(relres),
        relres = mean(relres)
      )
  }

  list(genes = gene_qc, bulks = bulk_qc)

}
