### load/define required inputs/outputs ###

# load input data 
out_fast_csre_noboot <- readRDS("../testdata/out_fast_csre_noboot.rds")
out_csre_noboot <- readRDS("../testdata/out_csre_noboot.rds")
out_csre <- readRDS("../testdata/out_csre.rds")
out_fit_tissue_noboot <- readRDS("../testdata/out_fit_tissue_noboot.rds")
out_fit_tissue <- readRDS("../testdata/out_fit_tissue.rds")

# get single cell data
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

# get sc mapping
mapping <- tibble::tibble(
  cell_id = colnames(test_sclibrary),
  celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
)

###


### TESTS ###

test_that("fast_specific_expression_regulation produces 
  correct output", {

  test_fast_csre_noboot <- suppressMessages(
    fast_specific_expression_regulation(
      out_fit_tissue_noboot,
      test_sclibrary,
      mapping,
      show_pb = FALSE
    )
  )

  expect_equal(test_fast_csre_noboot, out_fast_csre_noboot)

})

test_that("specific_expression_regulation produces 
  correct output without bootstrapping", {

  test_csre_noboot <- suppressMessages(
    specific_expression_regulation(
      out_fit_tissue_noboot,
      test_sclibrary,
      mapping
    )
  )

  expect_equal(test_csre_noboot, out_csre_noboot)

})

test_that("specific_expression_regulation produces 
  correct output with bootstrapping", {

  test_csre <- suppressMessages(
    specific_expression_regulation(
      out_fit_tissue,
      test_sclibrary,
      mapping
    )
  )

  expect_equal(test_csre, out_csre)

})

###