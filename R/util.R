## namespace: start
#' @importFrom dplyr %>%
#' @importFrom rlang := .data
## namespace: end
NULL

#' Filter two matrices to their common genes
#'
#' @description
#' Restricts rows of two matrices to the set of common genes and brings
#' them in the same ordering. Additionally it will report if less than
#' 20% of the genes are in common.
#'
#' @param mat_A Matrix A.
#' @param mat_B Matrix B.
#'
#' @return A list containing the input A and B but restricted to the set of
#'    common genes.
#'
#' @examples
#' set.seed(42)
#' ngenes <- 4
#' nbulks <- 3
#' ncells <- 5
#' bulks <- matrix(runif(ngenes * nbulks), nrow = ngenes, ncol = nbulks)
#' sc <- matrix(runif((ngenes + 1) * ncells), nrow = ngenes + 1, ncol = ncells)
#' bulknames <- paste0(rep("bulk_", nbulks), seq.int(1, nbulks))
#' cellnames <- paste0(rep("cell_", ncells), seq.int(1, ncells))
#' genenames_bulks <- paste0(rep("gene_", ngenes), seq.int(1, ngenes))
#' genenames_sc <- paste0(rep("gene_", ngenes + 1), seq.int(1, ngenes + 1))
#' rownames(bulks) <- genenames_bulks
#' rownames(sc) <- genenames_sc
#' colnames(bulks) <- bulknames
#' colnames(sc) <- cellnames
#'
#' common <- tissueResolver:::common_genes(bulks, sc)
#' bulks_common <- common$A
#' sc_common <- common$B
#'
common_genes <- function(
  mat_A,
  mat_B
) {

  if (is.null(rownames(mat_A)) || is.null(rownames(mat_B))) {
    stop("common_genes: matrices A and B must have row names (= gene names)")
  }
  genelist <- intersect(rownames(mat_A), rownames(mat_B))
  if (length(genelist) < 1) {
    stop("common_genes: no genes in common")
  } else if (
    length(genelist) < 0.2 * nrow(mat_A) || length(genelist) < 0.2 * nrow(mat_B)
  ) {
    warning("common_genes: less than 20% common genes!")
  }
  list(
    A = as.matrix(mat_A[genelist, ]),
    B = as.matrix(mat_B[genelist, ])
  )

}

#' Promote a vector to a matrix
#'
#' @description
#' Promotes a vector to a matrix. If the input is already a matrix, this matrix
#' is returned.
#'
#' @param vec A matrix or vector.
#'
#' @return Promoted matrix.
#'
promote_matrix <- function(
  vec
) {

  if (is.vector(vec)) {
    as.matrix(vec)
  } else if (is.matrix(vec)) {
    vec
  } else {
    stop("promote_matrix: input is neither vector nor matrix")
  }

}

#' Check input dimensions of bulk and single cell data
#'
#' @description
#' Checks if the number of genes of bulk data and the single cell library match.
#'
#' @param x A matrix representing the collection of bulk samples.
#' @param y A matrix representing the single cell library.
#'
#' @examples
#' set.seed(42)
#' ngenes <- 4
#' nbulks <- 3
#' ncells <- 5
#' bulks <- matrix(runif(ngenes * nbulks), nrow = ngenes, ncol = nbulks)
#' sc <- matrix(runif((ngenes + 1) * ncells), nrow = ngenes + 1, ncol = ncells)
#' bulknames <- paste0(rep("bulk_", nbulks), seq.int(1, nbulks))
#' cellnames <- paste0(rep("cell_", ncells), seq.int(1, ncells))
#' genenames_bulks <- paste0(rep("gene_", ngenes), seq.int(1, ngenes))
#' genenames_sc <- paste0(rep("gene_", ngenes + 1), seq.int(1, ngenes + 1))
#' rownames(bulks) <- genenames_bulks
#' rownames(sc) <- genenames_sc
#' colnames(bulks) <- bulknames
#' colnames(sc) <- cellnames
#'
#' common <- tissueResolver:::common_genes(bulks, sc)
#' bulks_common <- common$A
#' sc_common <- common$B
#' tissueResolver:::check_input_dimensions(bulks_common, sc_common)
#'
check_input_dimensions <- function(
  y, x
) {

  if (
    (!is.matrix(x)) || (!is.matrix(y)) || (!is.numeric(x)) || (!is.numeric(y))
  ) {
    stop("bulk and sclibrary must both be numeric matrices")
  }
  if (nrow(y) != nrow(x)) {
    stop("incompatible dimensions of bulk and scdata")
  }

}

#' Check input for `fast_specific_expression_regulation()`
#'
#' @description
#' Checks if the input for computing specific expression is valid (accessed by
#' `fast_specific_expression_regulation()`)
#'
#' @param celldf A dataframe composed of tissuemodel returned by `fit_tissue()`
#'    and a mapping prescribing the mapping of cell_ids to cell types.
#' @param sclibrary A single cell data matrix, columns are individual cells and
#'    must be indexable by the `cell_id` contained in the `celldf`.
#' @param by The name of the categorical variable defining groups of cells
#' @param weight The name of the column that holds the fitted cell weights
#' @param bulk_id The name of the column that identifies the bulks
#' @param cell_id The name of the column that contains the ids for the cells
#'    that match the column names of `sclibrary`
#'
#' @return Nothing.
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
#' fitted_tissue <- fit_tissue(bulks,sc)
#'
#' mapping <- tibble(
#'   cell_id = cellnames,
#'   celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
#' )
#'
#' celldf <- fitted_tissue$tissuemodel %>% inner_join(mapping, by = "cell_id")
#'
#' tissueResolver:::check_input_se(
#'   celldf = celldf,
#'   sclibrary = sc,
#'   by = "celltype",
#'   weight = "weight",
#'   bulk_id = "bulk_id",
#'   cell_id = "cell_id"
#' )
#'
check_input_se <- function(
  celldf,
  sclibrary,
  by,
  weight,
  bulk_id,
  cell_id
) {

  if (!tibble::is_tibble(celldf)) {
    stop("celldf must be a tibble")
  }
  if (
    !all(sapply(
      celldf %>% dplyr::pull(cell_id),
      function(x) {
        x %in% colnames(sclibrary)
      }
    ))
  ) {
    stop("celldf contains weights of cells that are not in the count matrix.")
  }
  # in explained_expression by will be set to false
  # because we do not want it to be checked here
  if ((by != FALSE) && !by %in% names(celldf)) {
    stop(paste0("by-selection (\"", by, "\") is not a name of celldf"))
  }
  if (!weight %in% names(celldf)) {
    stop(paste0("weight (\"", weight, "\") is not a name of celldf"))
  }
  if (!bulk_id %in% names(celldf)) {
    stop(paste0("bulk_id (\"", bulk_id, "\") is not a name of celldf"))
  }
  if (!cell_id %in% names(celldf)) {
    stop(paste0("cell_id (\"", cell_id, "\") is not a name of celldf"))
  }

}

#' Convert a tibble or dataframe to a matrix
#'
#' @description
#' Convert a tibble or dataframe to a matrix to allows to call functions
#' expecting matrices to be called with tibbles.
#'
#' @param df A dataframe to be converted into matrix.
#' @param rows The name in the input dataframe that corresponds to the rows of
#'    the resulting matrix, defaults to `"gene"`.
#' @param cols The name in the input dataframe that corresponds to the columns
#'    of the resulting matrix, defaults to `"cell_id"`.
#' @param values The name in the input dataframe that will hold the entries of
#'    the resulting matrix, defaults to `"expression"`.
#'
#' @return A matrix storing the values of `df`.
#'
#' @examples
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
#'
#' fitted_tissue <- fit_tissue(bulks, sc)
#'
#' sc_df <- tissueResolver:::matrix_to_df(
#'   sc,
#'   rows = "gene",
#'   cols = "cell_id",
#'   values = "expression"
#' )
#'
#' sc_mat <- tissueResolver:::df_to_matrix(
#'   sc_df,
#'   rows = "gene",
#'   cols = "cell_id",
#'   values = "expression"
#' )
#'
df_to_matrix <- function(
  df,
  rows = "gene",
  cols = "cell_id",
  values = "expression"
) {
  # give each cell a seperate column storing the gene expression
  tmp <- df %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(cols),
      values_from = dplyr::all_of(values)
    )
  # discard the rows column and store as matrix
  r <- tmp %>% dplyr::select(-dplyr::all_of(rows)) %>% as.matrix()
  # give the rows of the matrix the names in rows
  rownames(r) <- tmp %>% dplyr::pull(dplyr::all_of(rows))
  # give the cols of the matrix the names in cols
  colnames(r) <- tmp %>% dplyr::select(-dplyr::all_of(rows)) %>% names()
  return(r)
}


#' Convert a matrix to a tibble
#'
#' @description
#' Convert a matrix to a tibble to allows to call functions expecting tibbles
#' to be called with dataframes.
#'
#' @param genematrix A matrix to convert.
#' @param rows The dedicated name in the resulting dataframe storing the row
#'    names of the input matrix, defaults to `"gene"`.
#' @param cols The dedicated name in the resulting dataframe storing the column
#'    names of the input matrix, defaults to `"cell_id"`.
#' @param values The dedicated name in the resulting dataframe that holds the
#'    entries of the input matrix, defaults to `"expression"`.
#'
#' @return A tibble storing the matrix values and names.
#'
#' @examples
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
#'
#' fitted_tissue <- fit_tissue(bulks, sc)
#'
#' sc_df <- tissueResolver:::matrix_to_df(
#'    sc,
#'    rows = "gene",
#'    cols = "cell_id",
#'    values = "expression"
#'   )
#'
#' sc_mat <- tissueResolver:::df_to_matrix(
#'    sc_df,
#'    rows = "gene",
#'    cols = "cell_id",
#'    values = "expression"
#' )
#'
matrix_to_df <- function(
  genematrix,
  rows = "gene",
  cols = "cell_id",
  values = "expression"
) {

  # create a tibble from the input matrix and
  # store the rownames in a seperate column
  res <- cbind(
    dplyr::tibble(gene_tmp = rownames(genematrix)),
    dplyr::as_tibble(genematrix)
    ) %>%
    # now store all the cells in one column
    tidyr::pivot_longer(-gene_tmp, names_to = cols, values_to = values)
  # assign the first column (storing the rownames of
  # the input matrix) the desired name given by rows
  names(res)[1] <- rows
  return(res)

}

#' Check input mapping for a `tissuemodel`
#'
#' @description
#' Checks if the mapping dataframe is a valid mapping for the cells in
#' `tissuemodel`.
#'
#' @param tissuemodel A `tissuemodel` as returned by `fit_tissue()`.
#' @param mapping A dataframe that maps cell_ids to celltypes.
#' @param by The name of the column that categorizes the mapping dataframe.
#' @param cell_id The name of the column in the mapping dataframe that names
#'    the single cells.
#'
#' @return Nothing.
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
#'
#' fitted_tissue <- tissueResolver:::fit_tissue_noboot(bulks, sc)
#'
#' mapping <- tibble(
#'   cell_id = cellnames,
#'   celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
#' )
#'
#' tissueResolver:::check_input_mapping(
#'   fitted_tissue, mapping, by = "celltype", cell_id = "cell_id"
#' )
#'
check_input_mapping <- function(
  tissuemodel,
  mapping,
  by,
  cell_id
) {

  if (!by %in% names(mapping)) {
    stop(paste0(
      "categorizing column by (=", by,
      ") is not a column name of the mapping dataframe"
    ))
  }
  if (
    !(cell_id %in% names(tissuemodel$tissuemodel)) ||
      !(cell_id %in% names(mapping))
  ) {
    stop(paste0(
      "cell_id column (", cell_id,
      ") is not a column name of both tissuemodel and the mapping dataframe"
    ))
  }
  if (
    length(intersect(
      tissuemodel$tissuemodel[[cell_id]],
      mapping[[cell_id]]
    )) == 0
  ) {
    stop("there are no common cell_ids in tissuemodel and mapping")
  }

}

#' Check if an object is a tissuemodel
#'
#' @description
#' Checks if the object is a tissuemodel as returned by `fit_tissue()`.
#'
#' @param tissuemodel object returned by `fit_tissue_noboot()` to be checked
#'
#' @return `TRUE` if it is a tissuemodel and `FALSE` if not.
#'
#' @examples
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
#'
#' fitted_tissue <- fit_tissue(bulks, sc)
#'
#' tissueResolver:::is_tissuemodel(fitted_tissue)
#'
is_tissuemodel <- function(
  tissuemodel
) {

  has_fields <- "tissuemodel" %in% names(tissuemodel) &&
    "fitted_genes" %in% names(tissuemodel)
  if (!has_fields) {
    return(FALSE)
  }
  length(tissuemodel$fitted_genes) > 0

}

#' Set up a progress bar
#'
#' @description
#' Initializes a progress bar and returns a function to update it.
#'
#' @param total The total number of steps for the progress bar.
#' @param min The starting value for the progress bar (defaults to `0`).
#' @param style The style of the progress bar (defaults to `3`).
#' @param ... Additional arguments passed to `utils::txtProgressBar()`.
#'
#' @return A function to update the progress bar.
#'
progress_bar <- function(
  total,
  min = 0,
  style = 3,
  ...
) {

  pb <- utils::txtProgressBar(min = min, max = total, style = style, ...)
  counter <- 0
  update <- function(n = 1) {
    counter <<- counter + n
    utils::setTxtProgressBar(pb, counter)
  }
  update

}