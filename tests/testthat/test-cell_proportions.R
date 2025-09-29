### load/define required inputs/outputs ###

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

# load reference inputs
out_fit_tissue_noboot <- readRDS("../testdata/out_fit_tissue_noboot.rds")
out_fit_tissue <- readRDS("../testdata/out_fit_tissue.rds")

# load reference outputs
out_prop_single <- readRDS("../testdata/out_prop_single.rds")
out_prop_noboot <- readRDS("../testdata/out_prop_noboot.rds")
out_prop <- readRDS("../testdata/out_prop.rds")

###



### TESTS ###

# check complete output of cell_proportions function without bootstrapping
test_that("cell_proportions_singlemodel produces correct output", {

  test_prop_single <- suppressMessages(
    cell_proportions_singlemodel(
      tissuemodel = out_fit_tissue_noboot,
      mapping = mapping
    )
  )
  expect_equal(test_prop_single, out_prop_single)

})

# check complete output of cell_proportions function without bootstrapping
test_that("cell_proportions without bootstrapping produces correct output", {

  test_prop_noboot <- suppressMessages(
    cell_proportions(
      tissuemodel = out_fit_tissue_noboot,
      mapping = mapping
    )
  )
  expect_equal(test_prop_noboot, out_prop_noboot)

})

# check complete output of cell_proportions function with bootstrapping
test_that("cell_proportions with bootstrapping produces correct output", {

  test_prop <- suppressMessages(
    cell_proportions(
      tissuemodel = out_fit_tissue,
      mapping = mapping
    )
  )
  expect_equal(test_prop, out_prop)

})

###