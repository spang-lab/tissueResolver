#' fits a given single cell library to a collection of bulk samples by
#' weighting the single cell gene expressions 
#' (the columns of the library matrix) in order to
#' optimally match the bulk gene expressions
#'
#'
#' @param bulkdata dataframe containing the gene expression of
#' a bulk sample in each column. Rownames are expected to be gene identifieres.
#' @param sclibrary dataframe containing a single cell profile in
#' each column. Rownames are expected to be gene identifieres
#'
#' @return a list containing
#' tissuemodel: contains the computed cell weights for each bulk sample
#' fitted_genes: the gene names of fitted genes
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
#' fitted_tissue <- tissueResolver:::fit_tissue_noboot(bulks,sc)
#'
fit_tissue_noboot <- function(bulkdata, sclibrary) {
  message(paste0("fitting ", ncol(bulkdata), " bulks..."))
  if (is_tibble(sclibrary)) {
    sclibrary <- df_to_matrix(sclibrary)
  }
  bulkdata <- promote_matrix(bulkdata)

  # take only common gene rows of bulk and single cell data
  minput <- common_genes(bulkdata, sclibrary)

  # y stores bulk data
  # x stores sc data
  y <- minput$A
  x <- minput$B
  check_input_dimensions(y, x)
  genes <- rownames(y)

  ncells <- ncol(x)

  # if the sc library/bulks have no column names,
  # use the column indices as cell/bulk ids instead
  if (is.null(colnames(x))) {
    cell_ids <- seq(1, ncol(x))
  } else {
    cell_ids <- colnames(x)
  }
  if (is.null(colnames(y))) {
    bulk_ids <- seq(1, ncol(y))
  } else {
    bulk_ids <- colnames(y)
  }

  # the single cell library x is weighted by the vector 
  # cell_weights in order to fit the i-th bulk sample 
  # vector y[,i]
  fn <- function(cell_weights, i) {
    r <- y[,i] - x %*% cell_weights
    return (sum(r^2))
  }

  # grad_fn returns the gradient of fn wrt. to the cell_weights
  # this is a vector of lenght of the number of common genes
  grad_fn <- function(cell_weights, i) {
    r <- y[,i] - x %*% cell_weights
    dr <- 2*r
    return(-t(dr) %*% x)

  }

  result <- tibble()

  # print status depending on bulk samples
  pb <- progress::progress_bar$new(
                       format = "[:bar] :percent fitting :what",
                       total = ncol(y)
                     )
  for (ibulk in 1:ncol(y)) {
    pb$tick(tokens = list(what = bulk_ids[ibulk]))

    # perform the actual fitting by calculating the optimal cell_weights,
    # i.e., minimize the squared l2 loss
    cell_weights <- rep(1.0, ncells)

    # note the box constraint, demanding the cell weights to be non-negative
    res <- optim(
              cell_weights,
              fn,
              grad_fn,
              lower = rep(0, ncells),
              method = "L-BFGS-B",
              ibulk
            )
    # message status concerning convergence of optimizer
    convergence <- res$convergence
    if(convergence == 1){
      stop("iteration limit has been reached")
    }
    if(convergence == 51){
      stop(res$message)
    }
    if(convergence == 52){
      stop(res$message)
    }

    # retrieve the optimal parameters
    r <- res$par

    # for each bulk, result stores the so called tissuemodel, i.e.,
    # the cells and corresponding weights whenever strictly positive
    result <- rbind(result,
                    tibble(
                      bulk_id = bulk_ids[ibulk],
                      cell_id = cell_ids[r > 0.0],
                      weight  = r[r > 0.0]))
    message("done.")
  }


  return (list(
    tissuemodel = result,
    fitted_genes = genes
    ))
}


#' Fits a given single cell library to a collection of bulk samples by
#' weighting the single cell gene expressions
#' (the columns of the library matrix) in order to
#' optimally match the bulk gene expressions.
#' Allows for bootstrapping specified percentage of cells.
#'
#' @param bulkdata dataframe containing the gene expression of
#' a bulk sample in each column. Rownames are expected to be gene identifieres
#' @param sclibrary dataframe containing a single cell profile in
#' each column. Rownames are expected to be gene identifieres
#' @param bootstrap if set to FALSE fit tissue without bootstrapping.
#' If set to true select bootstrap_pctcells perecent of cells
#' and perform bootstrap_nruns of sampling cells and fitting these to bulks
#' @param bootstrap_nruns number of bootstrap runs
#' @param bootstrap_pctcells integer in [0,100] which prescribes the percentage of cells
#' sampled from the single cell library in each bootstrap run
#'
#' @return if bootstrap is FALSE, returns a list with entries
#' - tissuemodel: a data frame which holds all weights for each ' selected cell (id) and the bulk id that it was fitted to.
#' - fitted_genes: a list of fitted gene names
#' - bootstrap: FALSE
#' If bootstrap is TRUE it returns a list containing
#' - tissuemodels: a list of tissuemodels of every bootstrap run
#' - fitted_genes: a list of fitted gene names
#' - result$tissuemodels[[irun]] <- this_model$tissuemodel
#' - bootstrap: TRUE
#' - bootstrap_nruns: number of bootstrap samples
#' - bootstrap_pctcells: percentage of all single cells that were used in every bootstrap run
#' 
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
#' fitted_tissue <- fit_tissue(bulks,sc, bootstrap = TRUE, bootstrap_nruns = 5, bootstrap_pctcells = 50)
#' @export
#'
#
fit_tissue <- function(
                bulkdata,
                sclibrary,
                bootstrap = FALSE,
                bootstrap_nruns = 50,
                bootstrap_pctcells = 10) {
  if (bootstrap == TRUE){
    if( ! is.numeric(bootstrap_nruns) ){
      stop("bootstrap_nruns is not a number!")
    }
    if( ! is.numeric(bootstrap_pctcells) ){
      stop("bootstrap_pctcells is not a number!")
    }
    # select only a certain percentage of cells for bootstrapping
    ncells <- round(ncol(sclibrary) * bootstrap_pctcells / 100)
    result <- list(bootstrap = TRUE,
                   bootstrap_nruns = bootstrap_nruns,
                   bootstrap_pctcells = bootstrap_pctcells,
                   bootstrap_ncells = ncells,
                   tissuemodels = list(),
                   fitted_genes = c()
                   )
    for( irun in 1:bootstrap_nruns ){
      message(paste0("bootstrap run ", irun, " of ", bootstrap_nruns, " (using ", ncells, " cells)"))
      # sample ncells from the sc library (sample is different in each run)
      take_cells <- sample(colnames(sclibrary), ncells, replace=FALSE)

      # fit this sampled library to bulk
      this_model <- fit_tissue_noboot(bulkdata, sclibrary[,take_cells])

      # for each run store fitted_genes and tissuemodels
      result$fitted_genes <- this_model$fitted_genes
      result$tissuemodels[[irun]] <- this_model$tissuemodel
      result$bootstrap <- TRUE
      result$bootstrap_nruns <- bootstrap_nruns
      result$bootstrap_pctcells <- bootstrap_pctcells
    }
    return(result)
  } else {
    tm <- fit_tissue_noboot(bulkdata, sclibrary)
    tm$bootstrap <- FALSE
    return(tm)
  }
}
