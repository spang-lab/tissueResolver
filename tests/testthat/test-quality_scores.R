### load/define required inputs/outputs ###

# bulk data
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

# load reference inputs
out_csre_noboot <- readRDS("../testdata/out_csre_noboot.rds")
out_csre <- readRDS("../testdata/out_csre.rds")

# load reference outputs
out_qc_noboot <- readRDS("../testdata/out_qc_noboot.rds")
out_qc <- readRDS("../testdata/out_qc.rds")

###



### TESTS ###

test_that("quality_scores produces correct output without bootstrapping", {

  expect_warning(
    test_qc_noboot <- quality_scores(
      out_csre_noboot,
      test_bulkdata
      )
  )

  expect_equal(test_qc_noboot, out_qc_noboot)

})

test_that("quality_scores produces correct output with bootstrapping", {

  expect_warning(
    test_qc <- quality_scores(
      out_csre,
      test_bulkdata
      )
  )

  expect_equal(test_qc, out_qc)

})

###