### load/define required inputs/outputs ###

# load inputs
tm_noboot <- readRDS("../testdata/out_fit_tissue_noboot.rds")

###



### TESTS ###

# common_genes function
test_that("common_genes works as expected", {

  # define test matrices
  matA <- matrix(
    data = 1:6,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("sample", 1:2)
    )
  )

  matB <- matrix(
    data = 7:12,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 2:4),
      paste0("cell", 1:2)
    )
  )

  # expected output
  expected_out <- list(
    A = matA[c("gene2", "gene3"), ],
    B = matB[c("gene2", "gene3"), ]
    )

  # test function
  expect_equal(common_genes(matA, matB), expected_out)

})

test_that("common_genes handles no common genes", {

  # define test matrices
  matA <- matrix(
    data = 1:6,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("sample", 1:2)
    )
  )

  matB <- matrix(
    data = 7:12,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 4:6),
      paste0("cell", 1:2)
    )
  )

  # test function
  expect_error(
    common_genes(matA, matB),
    "common_genes: no genes in common"
    )

})

test_that("common_genes indicates low number of common genes", {

  # define test matrices
  matA <- matrix(
    data = 1:20,
    nrow = 10,
    ncol = 2,
    dimnames = list(
      paste0("gene", 1:10),
      paste0("sample", 1:2)
    )
  )

  matB <- matrix(
    data = 21:40,
    nrow = 10,
    ncol = 2,
    dimnames = list(
      paste0("gene", 10:19),
      paste0("cell", 1:2)
    )
  )

  # test function
  expect_warning(
    common_genes(matA, matB),
    "common_genes: less than 20% common genes!"
    )

})

# common_genes function
test_that("common_genes always outputs a list of matrices", {

  # define test matrices
  matA <- matrix(
    data = 1:6,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("sample", 1:2)
    )
  )

  matB <- matrix(
    data = 7:12,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 2:4),
      paste0("cell", 1:2)
    )
  )

  # transform to data.frame
  dfA <- as.data.frame(matA)
  dfB <- as.data.frame(matB)

  # test for list of 2 matrices
  for (input_pair in list(
    list(matA, matB),
    list(dfA, matB),
    list(matA, dfB),
    list(dfA, dfB)
  )) {
    out <- common_genes(input_pair[[1]], input_pair[[2]])
    expect_type(out, "list")
    expect_length(out, 2)
    expect_true(is.matrix(out[[1]]))
    expect_true(is.matrix(out[[2]]))
  }

  # test if outputs are equal
  expect_equal(
    common_genes(matA, matB),
    common_genes(dfA, dfB)
  )

})
#

# promote_matirx #
test_that("promote_matrix works as expected", {

  # define input vector
  in_vec <- 1:3

  # define expected output
  out_mat <- matrix(1:3, nrow = 3, ncol = 1)

  # as data.frame
  in_df <- as.data.frame(in_vec)

  # test function
  expect_equal(promote_matrix(in_vec), out_mat)
  expect_equal(promote_matrix(out_mat), out_mat)
  expect_true(is.matrix(promote_matrix(in_vec)))
  expect_true(is.matrix(promote_matrix(out_mat)))
  expect_error(
    promote_matrix(in_df),
    "promote_matrix: input is neither vector nor matrix"
    )

})
#

# check_input_dimensions #
test_that("check_input_dimensions works as expected", {

  # define test matrices
  matA <- matrix(
    data = 1:6,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("sample", 1:2)
    )
  )

  matB <- matrix(
    data = 7:12,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 2:4),
      paste0("cell", 1:2)
    )
  )

  matC <- matrix(
    data = 13:20,
    nrow = 4,
    ncol = 2,
    dimnames = list(
      paste0("gene", 1:4),
      paste0("cell", 1:2)
    )
  )

  matD <- matrix(
    data = c("a", "b", "c", "d", "e", "f"),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      paste0("gene", 2:4),
      paste0("cell", 1:2)
    )
  )

  # correct inputs -> no error/message
  expect_silent(check_input_dimensions(matA, matB))
  expect_error(
    check_input_dimensions(matA, matC),
    "incompatible dimensions of bulk and scdata"
    )
  expect_error(
    check_input_dimensions(matA, matC),
    "incompatible dimensions of bulk and scdata"
    )

})
#

# check_input_se #
test_that("check_input_se works as expected", {

  # single cell inputs
  sclibrary <- matrix(
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

  # mapping
  mapping <- tibble::tibble(
    cell_id = colnames(sclibrary),
    celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
  )

  # add grps to tm and get actual models
  celldf <- tm_noboot$tissuemodel %>%
    dplyr::inner_join(mapping, by = "cell_id")

  # function variables that stay consistent
  inputs <- list(
    by = "celltype",
    weight = "weight",
    bulk_id = "bulk_id",
    cell_id = "cell_id"
  )
    
  # correct inputs -> no error/message
  expect_silent(do.call(check_input_se, c(list(celldf, sclibrary), inputs)))

  # wrong cell_id
  celldf_wrong1 <- celldf
  celldf_wrong1$cell_id[1] <- "wrong_id"
  expect_error(
    do.call(check_input_se, c(list(celldf_wrong1, sclibrary), inputs)),
    "celldf contains weights of cells that are not in the count matrix."
  )

  # not a tibble
  celldf_wrong2 <- as.data.frame(celldf)
  expect_error(
    do.call(check_input_se, c(list(celldf_wrong2, sclibrary), inputs)),
    "celldf must be a tibble"
  )

  # wrong/missing column names
  for(wrong_input in list(
    list(
      by = "wrong_col",
      weight = "weight",
      bulk_id = "bulk_id",
      cell_id = "cell_id"
      ),
    list(
      by = "celltype",
      weight = "wrong_col",
      bulk_id = "bulk_id",
      cell_id = "cell_id"
      ),
    list(
      by = "celltype",
      weight = "weight",
      bulk_id = "wrong_col",
      cell_id = "cell_id"
      ),
    list(
      by = "celltype",
      weight = "weight",
      bulk_id = "bulk_id",
      cell_id = "wrong_col"
      )
    )
  ) {
    expect_error(
      do.call(check_input_se, c(list(celldf, sclibrary), wrong_input)),
      '(by-selection|weight|bulk_id|cell_id) \\(".+"\\) is not a name of celldf'
    )
  }

})
#

# df_to_matrix #
test_that("df_to_matrix works as expected", {

  # define test tibble
  tbl <- tibble::tibble(
  gene = paste0("Gene", 1:5) %>% rep(each = 3),
  cell_id = paste0("Cell", 1:3) %>% rep(times = 5),
  expression = c(
    # c1 c2 c3
      0, 1, 0,  # Gene1
      1, 0, 0,  # Gene2
      1, 1, 0,  # Gene3
      0, 1, 1,  # Gene4
      0, 0, 1   # Gene5
    )
  )

  # convert to data.frame
  df <- as.data.frame(tbl)

  # define expected output
  out_mat <- matrix(
    c(
    # c1 c2 c3
      0, 1, 0,  # Gene1
      1, 0, 0,  # Gene2
      1, 1, 0,  # Gene3
      0, 1, 1,  # Gene4
      0, 0, 1   # Gene5
    ),
    nrow = 5,
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      paste0("Gene", 1:5),
      paste0("Cell", 1:3)
    )
  )

  # test function
  expect_equal(df_to_matrix(tbl), out_mat)
  expect_equal(df_to_matrix(df), out_mat)

})
#

test_that("df_to_matrix works as expected", {

  # input matrix
  mat <- matrix(
    c(
    # c1 c2 c3
      0, 1, 0,  # Gene1
      1, 0, 0,  # Gene2
      1, 1, 0,  # Gene3
      0, 1, 1,  # Gene4
      0, 0, 1   # Gene5
    ),
    nrow = 5,
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      paste0("Gene", 1:5),
      paste0("Cell", 1:3)
    )
  )

  # define expected output
  out_tbl <- tibble::tibble(
  gene = paste0("Gene", 1:5) %>% rep(each = 3),
  cell_id = paste0("Cell", 1:3) %>% rep(times = 5),
  expression = c(
    # c1 c2 c3
      0, 1, 0,  # Gene1
      1, 0, 0,  # Gene2
      1, 1, 0,  # Gene3
      0, 1, 1,  # Gene4
      0, 0, 1   # Gene5
    )
  )

  # test function
  expect_equal(matrix_to_df(mat), out_tbl)

})
#

# check_input_mapping #
test_that("check_input_mapping works as expected", {

  # mapping
  mapping <- tibble::tibble(
    cell_id = paste0("Cell", 1:5),
    celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
  )

  # wrong_mappings 
  mapping_wrong_col_by <- mapping
  names(mapping_wrong_col_by)[2] <- "wrong_col"

  mapping_wrong_col_id <- mapping
  names(mapping_wrong_col_id)[1] <- "wrong_col"

  mapping_wrong_ids <- mapping
  mapping_wrong_ids$cell_id[1:5] <- "wrong_id"

  # correct inputs -> no error/message
  expect_silent(
    check_input_mapping(tm_noboot, mapping, by = "celltype", cell_id = "cell_id")
  )

  # faulty mappings
  expect_error(
    check_input_mapping(tm_noboot, mapping_wrong_col_by, by = "celltype", cell_id = "cell_id"),
    "categorizing column by \\(=celltype\\) is not a column name of the mapping dataframe"
  )

  expect_error(
    check_input_mapping(tm_noboot, mapping_wrong_col_id, by = "celltype", cell_id = "cell_id"),
    "cell_id column \\(cell_id\\) is not a column name of both tissuemodel and the mapping dataframe"
  )

  expect_error(
    check_input_mapping(tm_noboot, mapping_wrong_ids, by = "celltype", cell_id = "cell_id"),
    "there are no common cell_ids in tissuemodel and mapping"
  )

})
#

# is tissuemodel #
test_that("is_tissuemodel works as expected", {

  # correct input
  expect_true(is_tissuemodel(tm_noboot))

  # correct but missing fitted_genes
  tm_wrong1 <- tm_noboot
  tm_wrong1$fitted_genes <- NULL
  expect_false(is_tissuemodel(tm_wrong1))

  # correct but missing tissuemodel
  tm_wrong2 <- tm_noboot
  tm_wrong2$tissuemodel <- NULL
  expect_false(is_tissuemodel(tm_wrong2))

  # correct but fitted_genes are length 0
  tm_wrong3 <- tm_noboot
  tm_wrong3$fitted_genes <- character(0)
  expect_false(is_tissuemodel(tm_wrong3))

  # wrong inputs
  expect_false(is_tissuemodel(list(a = 1, b = 2)))
  expect_false(is_tissuemodel(data.frame(a = 1:3, b = 4:6)))

})
#

###