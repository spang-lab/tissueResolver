### load/define required inputs/outputs ###

# load inputs
csre <- readRDS("../testdata/out_csre.rds")

# inputs
grouping <- tibble::tibble(
  "bulk_id" = c("Bulk1", "Bulk2", "Bulk3", "Bulk4"),
  "group" = c("group_1", "group_2", "group_1", "group_2")
)

# add group info to csre object
csre_grp <- csre %>%
  dplyr::inner_join(grouping, by = c("bulk_id"))

###



### TESTS ###

test_that("plot_csre output is consistent (example settings)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    heatmapviz = "relchange",
    barplotviz = "relative",
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    addexpressionlevel = TRUE
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (default settings)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB")
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_default",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (heatmapviz = logdiffreg)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    heatmapviz = "logdiffreg"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_logdiffreg",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (heatmapviz = relregchange)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    heatmapviz = "relregchange"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_relregchange",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (barplotviz = log)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    barplotviz = "log"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_log",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (barplotviz = logdiff)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    barplotviz = "logdiff"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_logdiff",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (barplotviz = relative)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    barplotviz = "relative"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_relative",
    test_csre_plot
  )
  sink()

})

test_that("plot_csre output is consistent (barplotviz = relativeB)", {
  
  # tests
  test_csre_plot <- plot_csre(
    csre = csre_grp,
    groupA = "group_1",
    groupB = "group_2",
    ctypes = c("CTA", "CTB"),
    barplotviz = "relativeB"
  )
  
  sink(tempfile())
  vdiffr::expect_doppelganger(
    "plot_csre_relativeB",
    test_csre_plot
  )
  sink()

})

###