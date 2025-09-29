### load/define required inputs/outputs ###

# load inputs
prop_noboot <- readRDS("../testdata/out_prop_noboot.rds")
prop <- readRDS("../testdata/out_prop.rds")

# define grouping
grouping <- tibble::tibble(
  "bulk_id" = c("Bulk1", "Bulk2", "Bulk3", "Bulk4"),
  "group" = c("group_1", "group_2", "group_1", "group_2")
)

###



### TESTS ###

test_that("prop_visualization works as expected without bootsrapping", {

  # tests
  test_prop_plot <- prop_visualization(
    props = prop_noboot,
    grouping = grouping,
    groupA = "group_1",
    groupB = "group_2"
  )

  vdiffr::expect_doppelganger(
    "prop_visualization_noboot",
    test_prop_plot
  )
})

test_that("prop_visualization works as expected", {

  # tests
  test_prop_plot <- prop_visualization(
    props = prop,
    grouping = grouping,
    groupA = "group_1",
    groupB = "group_2"
  )

  vdiffr::expect_doppelganger(
    "prop_visualization",
    test_prop_plot
  )
})

test_that("prop_visualization works as expected (added grouping)", {

  # add groups to props
  prop_grp <- prop %>%
    dplyr::inner_join(grouping, by = "bulk_id")

  # tests
  test_prop_plot <- prop_visualization(
    props = prop_grp,
    grouping = grouping,
    groupA = "group_1",
    groupB = "group_2"
  )

  vdiffr::expect_doppelganger(
    "prop_visualization",
    test_prop_plot
  )
})

###