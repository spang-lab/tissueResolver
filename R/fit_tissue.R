#' fits a given single cell library to a collection of bulk samples by
#' weighting the single cell gene expressions
#' (the columns of the library matrix) in order to
#' optimally match the bulk gene expressions
#' The bulks are fit in parallel on multiple cores
#'
#'
#' @param bulkdata matrix containing the gene expression of
#' a bulk sample in each column. Rownames are expected to be gene identifieres.
#' @param sclibrary tibble or matrix containing a single cell profile in
#' each column. Rownames are expected to be gene identifieres
#' @param maxit the iteration limit used in the optimization algorithm, defaults to 2e4. This is very dataset specific,
#' so if the provided value is not sufficient we advise to increase this value.
#' @param ncores the number of cores for parallel computation
#'
#' @return a list containing
#' -tissuemodel: contains the computed cell weights for each bulk sample
#' -fitted_genes: the gene names of fitted genes
#' -log_store: logging info for residual across iterations
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
#' fitted_tissue <- tissueResolver:::fit_tissue_noboot(bulks, sc, 2e3, 2)
#'
fit_tissue_noboot <- function(bulkdata, sclibrary, maxit = 2e3, ncores = 1) {
  message(paste0("fitting ", ncol(bulkdata), " bulks..."))
  if (tibble::is_tibble(sclibrary)) {
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
  fn <- function(cell_weights, bulk_id) {
    r <- y[, bulk_id] - x %*% cell_weights
    return(sum(r^2))
  }

  # grad_fn returns the gradient of fn wrt. to the cell_weights
  # this is a vector of lenght of the number of common genes
  grad_fn <- function(cell_weights, bulk_id) {
    r <- y[, bulk_id] - x %*% cell_weights
    dr <- 2 * r
    return(-t(dr) %*% x)
  }
  # function for fitting a single bulk, will be called in parallel below
  fit_one_bulk <- function(bulk_id) {

    # perform the actual fitting by calculating the optimal cell_weights,
    # i.e., minimize the squared l2 loss
    cell_weights <- rep(1.0, ncells)

    # note the box constraint, demanding the cell weights to be non-negative
    # log for iteration info on residuals and convergence
    out_log <- capture.output(res <- optim(
      cell_weights,
      fn,
      grad_fn,
      lower = rep(0, ncells),
      method = "L-BFGS-B",
      control = list(maxit = maxit, trace = 1),
      bulk_id
    ))
    # message status concerning convergence of optimizer
    convergence <- res$convergence
    if (convergence == 1) {
      stop(paste0("iteration limit of", maxit, "has been reached"))
    }
    if (convergence == 51) {
      stop(res$message)
    }
    if (convergence == 52) {
      stop(res$message)
    }

    # retrieve the optimal parameters
    r <- res$par

    # for each bulk, result stores the so called tissuemodel, i.e.,
    # the cells and corresponding weights whenever strictly positive
    result <- tibble(
      bulk_id = bulk_id,
      cell_id = cell_ids[r > 0.0],
      weight  = r[r > 0.0]
    )

    # for each bulk store residuals for all iterations
    log_store <- tibble(
      bulk_id = bulk_id,
      log = list(out_log)
    )

    message(paste0("Fitting bulk ", bulk_id, " done."))

    return(list(
      result = result,
      log_store = log_store
    ))
  }
  # parallel execution of fitting bulks on ncores
  # bind all results together, note that mclapply keeps order of bulks
  all_fits <- mclapply(bulk_ids, fit_one_bulk, mc.cores = min(ncores, length(bulk_ids)))

  result <- tibble()
  log_store <- tibble()

  for (i in 1:ncol(y)) {
    result <- rbind(
      result,
      all_fits[[i]]$result
    )

    log_store <- rbind(
      log_store,
      all_fits[[i]]$log_store
    )
  }

  return(list(
    tissuemodel = result,
    fitted_genes = genes,
    log_store = log_store
  ))
}


#' Fits a given single cell library to a collection of bulk samples by
#' weighting the single cell gene expressions
#' (the columns of the library matrix) in order to
#' optimally match the bulk gene expressions.
#' Allows for bootstrapping specified percentage of cells.
#' The bulks are fit in parallel on multiple cores
#'
#' @param bulkdata matrix containing the gene expression of
#' a bulk sample in each column. Rownames are expected to be gene identifieres
#' @param sclibrary tibble or matrix containing a single cell profile in
#' each column. Rownames are expected to be gene identifieres
#' @param maxit the iteration limit used in the optimization algorithm, defaults to 2e3
#' @param bootstrap if set to FALSE fit tissue without bootstrapping.
#' If set to true select bootstrap_pctcells perecent of cells
#' and perform bootstrap_nruns of sampling cells and fitting these to bulks
#' @param bootstrap_nruns number of bootstrap runs
#' @param bootstrap_pctcells integer in \[0,100\] which prescribes the percentage of cells
#' sampled from the single cell library in each bootstrap run
#' @param ncores the number of cores for parallel computation
#'
#' @return if bootstrap is FALSE, returns a list with entries
#' - tissuemodel: a data frame which holds all weights for each ' selected cell (id) and the bulk id that it was fitted to.
#' - fitted_genes: a list of fitted gene names
#' - bootstrap: FALSE
#' If bootstrap is TRUE it returns a list containing
#' - tissuemodels: a list of tissuemodels of every bootstrap run
#' - fitted_genes: a list of fitted gene names
#' - log_stores: a vector of log messages containing the residual per bulk and bootstrap run.
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
#' fitted_tissue <- fit_tissue(bulks, sc, bootstrap = TRUE, bootstrap_nruns = 5, bootstrap_pctcells = 50, ncores = 2)
#' @export
#'
fit_tissue <- function(bulkdata,
                           sclibrary,
                           maxit = 2e3,
                           bootstrap = FALSE,
                           bootstrap_nruns = 50,
                           bootstrap_pctcells = 10,
                           ncores = 1) {
  if (bootstrap == TRUE) {
    if (!is.numeric(bootstrap_nruns)) {
      stop("bootstrap_nruns is not a number!")
    }
    if (!is.numeric(bootstrap_pctcells)) {
      stop("bootstrap_pctcells is not a number!")
    }
    # select only a certain percentage of cells for bootstrapping
    ncells <- round(ncol(sclibrary) * bootstrap_pctcells / 100)
    result <- list(
      bootstrap = TRUE,
      bootstrap_nruns = bootstrap_nruns,
      bootstrap_pctcells = bootstrap_pctcells,
      bootstrap_ncells = ncells,
      tissuemodels = list(),
      log_stores = list(),
      fitted_genes = c()
    )
    for (irun in 1:bootstrap_nruns) {
      message(paste0("bootstrap run ", irun, " of ", bootstrap_nruns, " (using ", ncells, " cells)"))
      # sample ncells from the sc library (sample is different in each run)
      take_cells <- sample(colnames(sclibrary), ncells, replace = TRUE)

      # fit this sampled library to bulk
      this_model <- fit_tissue_noboot(bulkdata, sclibrary[, take_cells], maxit, ncores)

      # for each run store fitted_genes and tissuemodels
      result$fitted_genes <- this_model$fitted_genes
      result$tissuemodels[[irun]] <- this_model$tissuemodel
      result$log_stores[[irun]] <- this_model$log_store
      result$bootstrap <- TRUE
      result$bootstrap_nruns <- bootstrap_nruns
      result$bootstrap_pctcells <- bootstrap_pctcells
    }
    return(result)
  } else {
    tm <- fit_tissue_noboot(bulkdata, sclibrary, maxit, ncores)
    tm$bootstrap <- FALSE
    return(tm)
  }
}
