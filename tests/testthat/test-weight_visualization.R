### load required inputs/outputs ###

# load outputs
out_kde2d <- readRDS("../testdata/out_kde2d.rds")
out_den_dist <- readRDS("../testdata/out_den_dist.rds")

###



### TESTS ###

# kde2d_weighted #
  test_that("kde2d_weighted produces correct output", {
    # define inputs
    test_df <- data.frame(
      x = c(
        rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)
        ), 
      y = c(rep(x = 1:3, each = 3)), 
      c.type = c("A", "B", "C", "A", "C", "A", "B", "B", "C"),
      group_0 = c(1, 0.5, 0.25,  rep(0, 6)),
      group_1 = c(rep(0, 3), rep(1, 3), rep(0, 3)),
      group_2 = rev(c(rep(1, 3), rep(0, 6)))
    )

    test_kde2d <- kde2d_weighted(
      x = test_df$x,
      y = test_df$y,
      weight = test_df$group_0
    )

    expect_equal(test_kde2d, out_kde2d)

  })
#

# bulk_group_input #
  test_that("bulk_group_input works with named vector", {

    vec_bulkgrp <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_3")

    test_bulk_group_input <- bulk_group_input(
      bulkgroups = vec_bulkgrp
    )

    expect_equal(
      test_bulk_group_input,
      c("b1" = "group_1", "b2" = "group_2", "b3" = "group_3")
    )

  })

  test_that("bulk_group_input works with tibble", {

    tbl_bulkgrp <- tibble::tibble(
      "bulk_id" = c("b1", "b2", "b3"),
      "group" = c("group_1", "group_2", "group_3")
      )

    test_bulk_group_input <- bulk_group_input(
      bulkgroups = tbl_bulkgrp
    )

    expect_equal(
      test_bulk_group_input,
      c("b1" = "group_1", "b2" = "group_2", "b3" = "group_3")
    )

  })
#

# density distribution #
  test_that("density_distribution produces correct output", {

    # inputs
    sc_emb <- data.frame(
      "cell_id" = paste0("cell_", 1:9),
      "UMAP1" = rep(c(1, 2, 3), each = 3),
      "UMAP2" = rep(x = 1:3, each = 3),
      "celltype" =  c("A", "B", "C", "A", "C", "A", "B", "B", "C")
    )

    tm <- tibble::tibble(
      "fitted_genes" = NA,
      "tissuemodel" = tibble::tibble(
        "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
        "cell_id" = c(paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)),
        "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
      )
    )

    bulkgrp <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_1")

    # tests
    invisible(capture.output(
      test_den_dist <- density_distribution(
          scembedding = sc_emb,
          tissuemodel = tm$tissuemodel,
          bulkgroups = bulkgrp
        )
    ))
  
    expect_equal(test_den_dist, out_den_dist)

  })
#

# plot_element_embedding #
test_that("plot_element_embedding output is consistent", {

  # inputs
  sc_emb <- data.frame(
    "cell_id" = paste0("cell_", 1:9),
    "UMAP1" = rep(c(1, 2, 3), each = 3),
    "UMAP2" = rep(x = 1:3, each = 3),
    "celltype" =  c("A", "B", "C", "A", "C", "A", "B", "B", "C")
  )

  # tests
  test_embedding_plot <- plot_element_embedding(
    scembedding = sc_emb,
    colourby = "celltype"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_element_embedding",
    test_embedding_plot
  )
  sink()

})

test_that("plot_element_differential_density output is consistent", {

  # inputs
  sc_emb <- data.frame(
    "cell_id" = paste0("cell_", 1:9),
    "UMAP1" = rep(c(1, 2, 3), each = 3),
    "UMAP2" = rep(x = 1:3, each = 3),
    "celltype" =  c("A", "B", "C", "A", "C", "A", "B", "B", "C")
  )

  tm <- tibble::tibble(
    "fitted_genes" = NA,
    "tissuemodel" = tibble::tibble(
      "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
      "cell_id" = c(paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)),
      "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
    )
  )

  bulkgrp <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_1")

  invisible(capture.output(
    out_den_dist <- density_distribution(
      scembedding = sc_emb,
      tissuemodel = tm$tissuemodel,
      bulkgroups = bulkgrp
    )
  ))

  # tests
  test_diff_den_plot <- plot_element_differential_density(
    densities = out_den_dist,
    groupA = "group_1",
    groupB = "group_2"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_element_differential_density",
    test_diff_den_plot
  )
  sink()
  
})

test_that("plot_differential_density output is consistent", {

  # inputs
  sc_emb <- data.frame(
    "cell_id" = paste0("cell_", 1:9),
    "UMAP1" = rep(c(1, 2, 3), each = 3),
    "UMAP2" = rep(x = 1:3, each = 3),
    "celltype" =  c("A", "B", "C", "A", "C", "A", "B", "B", "C")
  )

  tm <- tibble::tibble(
    "fitted_genes" = NA,
    "tissuemodel" = tibble::tibble(
      "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
      "cell_id" = c(paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)),
      "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
    )
  )

  bulkgrp <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_1")

  invisible(capture.output(
    out_den_dist <- density_distribution(
      scembedding = sc_emb,
      tissuemodel = tm$tissuemodel,
      bulkgroups = bulkgrp
    )
  ))

  # tests
  test_diff_den_plot <- plot_differential_densities(
    density = out_den_dist,
    scembedding = sc_emb,
    groupA = "group_1",
    groupB = "group_2",
    colourby = "celltype"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_differential_density",
    test_diff_den_plot
  )
  sink()
  
})
#

###