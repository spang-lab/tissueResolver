### define required functions ###

# original fit_tissue functon
original_fit_tissue <- function(bulkdata,
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

###



### load/define required inputs/outputs ###

# test input data
test_bulkdata <- matrix(
  c(
    5, 3, 0, 1,
    4, 0, 2, 1,
    0, 2, 3, 5,
    0, 1, 4, 4,
    1, 0, 5, 3
    ),
  nrow = 5,
  ncol = 4,
  dimnames = list(
    paste0("Gene", 1:5),
    paste0("Bulk", 1:4)
  )
)

test_sclibrary <- matrix(
  c(
    1, 0, 0, 0, 0,
    1, 1, 0, 0, 0,
    0, 1, 1, 0, 0,
    0, 0, 1, 1, 0,
    0, 0, 0, 1, 1
    ),
  nrow = 5,
  ncol = 5,
  dimnames = list(
    paste0("Gene", 1:5),
    paste0("Cell", 1:5)
  )
)

# load pregenerated reference outputs
out_fit_tissue_noboot <- readRDS("../testdata/out_fit_tissue_noboot.rds")
out_fit_tissue <- readRDS("../testdata/out_fit_tissue.rds")

###



### TESTS ###

# check complete output of fit_tissue_noboot function
test_that("fit_tissue_noboot produces correct output", {

  test_fit_tissue_noboot <- suppressMessages(
    fit_tissue_noboot(
      test_bulkdata,
      test_sclibrary
    )
  )

  expect_equal(test_fit_tissue_noboot, out_fit_tissue_noboot)
  
})

# test_that("fit_tissue produces correct output", {

#   set.seed(123)
#   test_fit_tissue <- fit_tissue(
#     test_bulkdata,
#     test_sclibrary,
#     bootstrap = TRUE,
#     bootstrap_nruns = 5,
#     bootstrap_pctcells = 50
#   )
#   expect_equal(test_fit_tissue, out_fit_tissue)
  
# })

test_that("fit_tissue and fit_tissue_noboot return equal outputs", {

  # basic function
  test_basic <- suppressMessages(
    fit_tissue_noboot(
      test_bulkdata,
      test_sclibrary
    )
  )
  # new function
  test_wrapper <- suppressMessages(
    fit_tissue(
      test_bulkdata,
      test_sclibrary,
      bootstrap = FALSE
    )
  )

  expect_equal(test_wrapper, test_basic)

})

test_that("fit_tissue produces correct output with bootstrapping", {

  # original function
  set.seed(123)
  original_fit_tissue <- suppressMessages(
    fit_tissue(
      test_bulkdata,
      test_sclibrary,
      bootstrap = TRUE,
      bootstrap_nruns = 5,
      bootstrap_pctcells = 50
    )
  )

  # new function
  set.seed(123)
  test_fit_tissue <- suppressMessages(
    fit_tissue(
      test_bulkdata,
      test_sclibrary,
      bootstrap = TRUE,
      bootstrap_nruns = 5,
      bootstrap_pctcells = 50
    )
  )

  expect_equal(test_fit_tissue, original_fit_tissue)

})

###