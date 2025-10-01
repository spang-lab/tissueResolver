#' Get weighted 2D Kernel density estimation
#'
#' @description
#' This function is an adaption from `MASS::kde2d()`, using a weighted approach
#' to estimate the density of points in a 2D space.
#'
#' @param x The numeric x coordinate of the embedding.
#' @param y The numeric y coordinate of the embedding.
#' @param weight The numeric weight for each data point.
#' @param n_grid The integer number of grid points used for plotting. Defaults
#'    to `25`.
#' @param lims The numeric vector with a length of 4 `c(lower(x), upper(x),
#'    lower(y), upper(y))`, range of the plot, allowing to plot a subframe of
#'    the data. Defaults to 'c(range(x), range(y))'.
#' @param bandwidth The bandwidth for the KDE, inducing the smoothing for the
#'    density estimate. A numeric vector with length of 2. Defaults to
#'    `c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))/4`.
#'
#' @return A dataframe with 3 entries: `x` and `y` grid coordinates in long
#'    format and `z` (= weighted density) in long format.
#'
#' @examples
#' rm(list = ls()); gc()
#' library(tissueResolver)
#' library(ggplot2)
#'
#' set.seed(42)
#' n <- 10
#' tmp.frame <- data.frame(
#'   x = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   y = c(rep(x = 1:3, each = 3)),
#'   c.type = sample(c("A", "B", "C"), 9, replace = TRUE),
#'   group_0 = c(1, 0.5, 0.25,  rep(0, 6)),
#'   group_1 = c(rep(0, 3), rep(1, 3), rep(0, 3)),
#'   group_2 = rev(c(rep(1, 3), rep(0, 6)))
#' )
#'
#' weighted.density <- tissueResolver:::kde2d_weighted(
#'   x = tmp.frame$x
#'   , y = tmp.frame$y
#'   , weight = tmp.frame$group_0
#' )
#' # Note that the corresponding weight matrix is given by replicating
#' # the t(weight) across all rows of the grid
#'
#' p <- ggplot(
#'   data = tmp.frame
#'   , mapping = aes(
#'     x = x
#'     , y = y
#'     )
#'   ) +
#'   geom_tile(
#'     data = weighted.density
#'     , mapping = aes(
#'       x = x
#'       , y = y
#'       , fill = z
#'     )
#'   ) +
#'   geom_point()
#' print(p)
#'
kde2d_weighted <- function(
  x,
  y,
  weight,
  bandwidth = NULL,
  n_grid = 25,
  lims = c(range(x), range(y))
) {

  # number of data points:
  n_x <- length(x)

  if (length(y) != n_x) {
    stop("In 'kde2d_weighted': 'length(x)' and 'length(y)' differ.")
  }
  if (length(weight) != n_x) {
    stop("In 'kde2d_weighted': 'length(x)' and 'length(weight)' differ.")
  }

  grid_x <- seq(from = lims[1], to = lims[2], length = n_grid)
  grid_y <- seq(from = lims[3], to = lims[4], length = n_grid)

  if (is.null(bandwidth)) {
    bandwidth <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y)) / 4
  }

  # distance of each data point to each grid point scaled by the bandwidth,
  # for computation of kernel
  distance_x <- outer(grid_x, x, "-") / bandwidth[1]
  distance_y <- outer(grid_y, y, "-") / bandwidth[2]

  # weight matrix
  left_mat <- matrix(
    rep(weight, n_grid),
    nrow = n_grid,
    ncol = n_x,
    byrow = TRUE
  )

  # matrices for the kernel
  middle_mat <- matrix(stats::dnorm(distance_x), n_grid, n_x)
  right_mat <- matrix(stats::dnorm(distance_y), n_grid, n_x)

  # weighted kde differs from classical kde by multiplying the kernel with
  # weight matrix and dividing by the cumulated weights. For classical kde,
  # see Section 5.6 in Venables, W. N. and Ripley, B. D. (2002) Modern
  # Applied Statistics with S. Fourth edition. Springer.
  z <- (left_mat * middle_mat) %*% t(right_mat) /
    (sum(weight) * bandwidth[1] * bandwidth[2])

  # the returned dataframe stores the grid coordinates and the weighted kde
  ret <- data.frame(
    x = rep(grid_x, n_grid),
    y = rep(grid_y, each = n_grid),
    z = as.vector(z)
  )
  return(ret)
}

#' Format bulk groups input
#'
#' @description
#' An input wrapper for `density_distribution()` to accept both named vectors as
#' well as dataframes containing bulk groups.
#'
#' @param bulkgroups A vector or dataframe holding the bulk-group association.
#' @param bulkidcol If bulkgroups is a tibble, the name of the column that holds
#'    the bulk ids.
#' @param groupcol If bulkgroups is a tibble, the name of the column that holds
#'    the group ids.
#'
#' @return A named vector of group ids with corresponding bulk ids as names.
#'
#' @examples
#' # Case 1
#' bulkgroups <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_3")
#'
#' bulk_group_input_1 <- tissueResolver:::bulk_group_input(
#'    bulkgroups,
#'    bulkidcol = "bulk_id",
#'    groupcol = "group"
#' )
#'
#' # Case 2
#' library(tibble)
#' bulkgroups <- tibble(
#'   "bulk_id" = c("b1", "b2", "b3"),
#'   "group" = c("group_1", "group_2", "group_3")
#' )
#'
#' bulk_group_input_2 <- tissueResolver:::bulk_group_input(
#'    bulkgroups,
#'    bulkidcol = "bulk_id",
#'    groupcol = "group"
#' )
#'
bulk_group_input <- function(
  bulkgroups,
  bulkidcol = "bulk_id",
  groupcol = "group"
) {

  if (tibble::is_tibble(bulkgroups)) {
    if (!bulkidcol %in% names(bulkgroups)) {
      stop(paste0(
        bulkidcol, " is not a column of the bulk group definition dataframe"
      ))
    }
    if (!groupcol %in% names(bulkgroups)) {
      stop(paste0(
        groupcol, " is not a column of the bulk group definition dataframe"
      ))
    }
    tmp <- bulkgroups[[groupcol]]
    names(tmp) <- bulkgroups[[bulkidcol]]
    return(tmp)
  } else if (is.vector(bulkgroups)) {
    if (is.null(names(bulkgroups))) {
      stop(
        "bulk group input vector is unnamed. Expect bulk ids as names of 
          the group vector"
      )
    }
    return(bulkgroups)
  } else {
    stop("invalid input for bulkgroups")
  }

}

#' Get density distributions of tissuemodels
#'
#' @description
#' Computes the weighted 2D densities for cell populations in a tissuemodel
#' derived from `fit_tissue()`. If groups are given, the group-specific
#' densities are computed.
#'
#' @param scembedding A dataframe that holds the 2D embeddings of the single
#'    cell data used to compute the tissue model.
#' @param x The name of the x coordinate of the 2D embedding. Defaults to
#'    `"UMAP1"`.
#' @param y The name of the y coordinate of the 2D embedding. Defaults to
#'    `"UMAP2"`.
#' @param scids The name of the column holding the cell identifiers in the
#'    embedding data frame. Defaults to `"cell_id"`.
#' @param tissuemodel A tissuemodel derived from `fit_tissue()`.
#' @param model_cellids The name of the cell_id column in the model data frame.
#'    Defaults to `"cell_id"`.
#' @param model_bulkids The name of the bulk id column in the model data frame.
#'    Defaults to `"bulk_id"`.
#' @param bulkgroups An optional mapping for the bulks to additional groups
#'    (e.g. treated/untreated, etc).
#' @param bulk_groupids The names of the column that holds the group id within
#'    bulkgroups if that is a dataframe. Defaults to `"group"`.
#' @param verbose Prints more output if `TRUE` (default).
#'
#' @return A dataframe with a kde density z-value for every x and y coordinate
#'    and an additional group column.
#'
#' @examples
#' library(tibble)
#' set.seed(42)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#'
#' fitted_tissue <- tibble(
#'   "fitted_genes" = NA,
#'   "tissuemodel" = tibble(
#'     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#'     "cell_id" = c(
#'       paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)
#'     ),
#'     "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
#'   )
#' )
#'
#' bulkgroups <- c("b1" = "bgroup_1", "b2" = "bgroup_2", "b3" = "bgroup_1")
#'
#' density <- density_distribution(
#'   scembedding = scembedding,
#'   x = "UMAP1", y = "UMAP2", scids = "colnames",
#'   fitted_tissue$tissuemodel,
#'   bulkgroups = bulkgroups
#' )
#'
#' @export
#'
density_distribution <- function(
  scembedding,
  x = "UMAP1",
  y = "UMAP2",
  scids = "cell_id",
  tissuemodel,
  model_cellids = "cell_id",
  model_bulkids = "bulk_id",
  bulkgroups = NULL,
  bulk_groupids = "group",
  verbose = TRUE
) {

  # can deal with "full" tissuemodel, i.e. output of fit_tissue -> filter
  # only tissuemodel tibble
  if (is_tissuemodel(tissuemodel)) {
    tissuemodel <- tissuemodel$tissuemodel
    # can also deal with only tissuemodel tibble, then no filtering necessary
  } else if (!tibble::is_tibble(tissuemodel)) {
    stop("tissuemodel is not a tissuemodel or dataframe")
  }

  if (!model_bulkids %in% colnames(tissuemodel)) {
    stop(paste0("no bulk id group \"", model_bulkids, " in tissuemodel."))
  }
  if (is.null(bulkgroups)) {
    # if no bulkgroups are given
    # introduce the group vector with just one generic group "all"
    bulks <- unique(tissuemodel[[model_bulkids]])
    bulkgroups <- rep("all", length(bulks))
    names(bulkgroups) <- bulks
  } else {
    # if bulkgroups are given
    # format the mapping into a named vector
    bulkgroups <- bulk_group_input(bulkgroups, model_bulkids, bulk_groupids)
  }

  # sanity checks --------------------------------------------------------------
  if (!scids %in% colnames(scembedding)) {
    stop(
      paste0(
        "In 'density_distribution': there is no '", scids, "' column, ",
        " in 'scembedding'. Therefore I can not map weights from ",
        "'tissuemodel' to scRNA-seq profiles"
      )
    )
  }

  if (!model_cellids %in% colnames(tissuemodel)) {
    stop(
      paste0(
        "In 'density_distribution': there is no '", model_cellids, "' column, ",
        " in 'tissuemodel'. Therefore I can not map ",
        "weights from 'tissuemodel' to scRNA-seq profiles"
      )
    )
  }
  all.bulks <- names(bulkgroups)
  if (any(duplicated(all.bulks))) {
    stop(
      paste0(
        "In 'density_distribution': 'names(bulkgroups)' must be unique."
      )
    )
  }
  # only consider unique bulk ids which both are given in the tissuemodel
  # and in the bulkgroups vector
  fitted.bulks <- unique(tissuemodel[[model_bulkids]])
  bulks_intersect <- intersect(all.bulks, fitted.bulks)

  if (length(bulks_intersect) < 2) {
    stop(
      paste0(
        "In 'density_distribution': Only ",
        length(bulks_intersect), " bulks detected. Provide at least 2."
      )
    )
  }

  # end of sanity checks -------------------------------------------------------

  if (verbose) {
    print(
      paste0(
        "--> Detected a bulk intersect of ",
        length(bulks_intersect), ". Start to average per group."
      )
    )
  }
  # ATTENTION: Theoretically there could be different bulk ids in the
  # bulkgroup mapping compared to the tissuemodel's bulk ids
  unique.groups <- unique(bulkgroups)
  if (verbose) {
    print(
      paste0(
        "--> Detected ", length(unique.groups), " groups."
      )
    )
  }

  densities <- tibble::tibble()

  for (groupid in unique.groups){
    # filter the bulk ids from the boulkgroup mapping belonging to one group
    idx_this_groups <- which(bulkgroups == groupid)
    bulks_this_group <- names(bulkgroups)[idx_this_groups]

    # filter the bulk ids from the tissuemodel belonging to one bulk group
    idx_fitted_tissue <- which(
      tissuemodel[[model_bulkids]] %in% bulks_this_group
    )
    # store the part of the tissuemodel belonging to this group
    fits_for_this_group <- tissuemodel[idx_fitted_tissue, ]

    print(groupid)

    # stores for each cell in this group the mean across all weights
    # belonging to this cell
    weights_for_this_group <- sapply(
      X = unique(fits_for_this_group[[model_cellids]]),
      # function that computes the mean over all weights
      # of a certain cell belonging to this group
      FUN = function(cell) {
        # retrieve the weight(s) for a certain cell within this group
        idx_for_this_cell <- which(fits_for_this_group[[model_cellids]] == cell)
        weights_for_this_cell <- as.numeric(
          fits_for_this_group[idx_for_this_cell, ][["weight"]]
        )
        mean(weights_for_this_cell)
      }
    )
    weights_for_all_cells <- rep(0, nrow(scembedding))
    names(weights_for_all_cells) <- scembedding[[scids]]
    weights_for_all_cells[names(weights_for_this_group)] <- 
      weights_for_this_group

    # store for each group the density matrix of the kde, where the weight
    # vector is given by the mean cell weight for each cell within this group
    densities <- rbind(
      densities,
      tibble::as_tibble(kde2d_weighted(
        x = scembedding[[x]],
        y = scembedding[[y]],
        weight = weights_for_all_cells
      )) %>%
      dplyr::mutate(group = groupid))
  }
  return(densities)

}

#' Plots the single cell embedding
#'
#' @description
#' Plots each UMAP point from the single cell embedding. These may be coloured
#' (e.g. by the associated cell type).
#'
#' @param scembedding A dataframe that holds 2D embeddings of the single cell
#'    data used to compute the tissuemodel.
#' @param x The column name of the x coordinate for the 2D embedding in
#'    `scembedding`.
#' @param y The column name of the y coordinate for the 2D embedding in
#'    `scembedding`.
#' @param colourby If not `NULL`, colour the single cells by this column of
#'    `scembedding` (e.g. celltype).
#'
#' @return A ggplot plot object that can further be customized, viewed or saved
#'    in a file. For more details see `ggplot2`.
#'
#' @examples
#' library(ggplot2)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#'
#' embedding_plot <- tissueResolver:::plot_element_embedding(
#'   scembedding,
#'   x = "UMAP1",
#'   y = "UMAP2",
#'   colourby = "c.type"
#' )
#'
#' print(ggplot() +
#'   embedding_plot +
#'   coord_fixed())
#'
plot_element_embedding <- function(
  scembedding,
  x = "UMAP1",
  y = "UMAP2",
  colourby = NULL
) {

  if (is.null(colourby)) {
    return(
      ggplot2::geom_point(
        data = scembedding,
        ggplot2::aes(
          x = .data[[x]],
          y = .data[[y]]
        ),
        alpha = 0.5
      )
    )
  } else {
    if (!colourby %in% names(scembedding)) {
      stop(paste0(colourby, " is not a column of scembedding."))
    }
    return(
      ggplot2::geom_point(
        data = scembedding,
        ggplot2::aes(
          x = .data[[x]],
          y = .data[[y]],
          colour = forcats::as_factor(.data[[colourby]])
        ),
        alpha = 0.5
      )
    )
  }

}

#' Plot the differential density between two groups as a tile plot
#'
#' @description
#' The kde density at each grid point of the single cell embedding for each
#' bulkgroup provided by `density_distribution()` is used to assess the
#' difference in densities between two groups. This differential density is
#' then visualized in the kde grid space as tile plot. If group labels are
#' omitted, they are guessed (if only two groups exist).
#'
#' @param densities A dataframe returned by `density_distribution()`.
#' @param groupA The name of the group with positive weights. These are colored
#'    red, if weights are larger than in `groupB`.
#' @param groupB The name of the group with negative weights. These are colored
#'    blue, if weights are larger than in `groupA`.
#' @param x The column name of the x coordinate of the kde grid.
#' @param y The column name of the y coordinate of the kde grid.
#' @param z The column name of the z coordinate for the kde grid.
#' @param group The column name of the densities dataframe storing the
#'    bulkgroups.
#'
#' @return A ggplot plot object that can be further customized, viewed or saved
#'    in a file. For more details see `ggplot2`.
#'
#' @examples
#' library(tibble)
#' library(ggplot2)
#' set.seed(42)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#' fitted_tissue <- tibble(
#'   "fitted_genes" = NA,
#'   "tissuemodel" = tibble(
#'     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#'     "cell_id" = c(
#'       paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)
#'     ),
#'     "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
#'   )
#' )
#' bulkgroups <- c("b1" = "bgroup_1", "b2" = "bgroup_2", "b3" = "bgroup_1")
#' density <- density_distribution(
#'   scembedding = scembedding,
#'   x = "UMAP1", y = "UMAP2", scids = "colnames",
#'   fitted_tissue$tissuemodel,
#'   bulkgroups = bulkgroups
#' )
#'
#' element_differential_density <-
#'   tissueResolver:::plot_element_differential_density(
#'     density,
#'     groupA = "bgroup_1",
#'     groupB = "bgroup_2",
#'     x = "x",
#'     y = "y",
#'     z = "z",
#'     group = "group"
#' )
#'
#' print(ggplot() +
#'   element_differential_density +
#'   coord_fixed())
#'
plot_element_differential_density <- function(
  densities,
  groupA = NA,
  groupB = NA,
  x = "x",
  y = "y",
  z = "z",
  group = "group"
) {

  groups <- unique(densities[[group]])
  ngroups <- length(groups)

  # if there are only two groups in the densities df and one of the input groups
  # is NA, specify the NA group as the remaining group
  if (ngroups == 2) {
    if (is.na(groupA) && (!is.na(groupB))) {
      groupA <- dplyr::setdiff(groups, groupB)
    } else if (is.na(groupB) && (!is.na(groupA))) {
      groupB <- dplyr::setdiff(groups, groupA)
    }
  } else if (is.na(groupA) || is.na(groupB)) {
    stop("no groups given. explicitly set groupA and groupB.")
  } else if (groupA == groupB) {
    stop("two differential groups must not be equal.")
  }
  if ((!groupA %in% groups) || (!groupB %in% groups)) {
    stop("groups for differential densities are not in the defined groups.")
  }
  # compute the differential densities for the two given groups
  diffdensity <- densities %>% 
    # filter only the rows of the df belonging to the two groups
    dplyr::filter(.data[[group]] %in% c(groupA, groupB)) %>%
    # give each group a separate column
    tidyr::pivot_wider(
      names_from = dplyr::all_of(group),
      values_from = dplyr::all_of(z),
      names_prefix = "group_"
    ) %>%
    # give both groups the names specified by the input groups
    dplyr::rename(groupA = dplyr::all_of(paste0("group_", groupA))) %>%
    dplyr::rename(groupB = dplyr::all_of(paste0("group_", groupB))) %>%
    # store the differential density between groupA and groupB in the column
    # diff
    dplyr::mutate(diff = groupA - groupB) %>%
    # delete the seperate densities and keep only the differential density
    dplyr::select(dplyr::all_of(c(x, y, "diff")))
  return(
    # represent the density as coloured tiles
    ggplot2::geom_tile(
      data = diffdensity,
      ggplot2::aes(
        x = .data[[x]],
        y = .data[[y]],
        fill = diff
      )
    )
  )

}

#' Plot differential densities between two groups
#'
#' @description
#' Visualizes groups for which a tissue was fitted and their differences. The
#' difference in density between groupA and groupB is computed and plotted. If
#' group labels are omitted, they are guessed (if only two groups exist).
#' In addition, the single cell embedding and a contour plot of the density for
#' each bulk group are plotted.
#'
#' @param densitydata A dataframe returned by `density_distribution()`.
#' @param scembedding A dataframe or matrix that holds 2D embeddings of the
#'    single cell data used to compute the tissuemodel. Defaults to `NULL`. If
#'    `scembedding` is provided the single cell embedding is plotted.
#' @param groupA The name of the group with positive weights. These are colored
#'    red, if weights are larger than in `groupB`.
#' @param groupB The name of the group with negative weights. These are colored
#'    blue, if weights are larger than in `groupA`.
#' @param embedding_x The column name of the x coordinate for the 2D embedding
#'    in `scembedding`. Defaults to `"UMAP1"`.
#' @param embedding_y The column name of the y coordinate for the 2D embedding
#'    in `scembedding`. Defaults to `"UMAP2"`.
#' @param colourby defaults to `NULL`. If column name is provided, colour the
#'    single cells by this column of `scembedding` (e.g., celltype).
#'
#' @return A ggplot plot object that can be further customized, viewed or saved
#'    in a file. For more details see `ggplot2`.
#'
#' @examples
#' library(tibble)
#' library(ggplot2)
#' set.seed(42)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#' fitted_tissue <- tibble(
#'   "fitted_genes" = NA,
#'   "tissuemodel" = tibble(
#'     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#'     "cell_id" = c(
#'       paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)
#'     ),
#'     "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
#'   )
#' )
#' bulkgroups <- c("b1" = "bgroup_1", "b2" = "bgroup_2", "b3" = "bgroup_1")
#' density <- density_distribution(
#'   scembedding = scembedding,
#'   x = "UMAP1", y = "UMAP2", scids = "colnames",
#'   fitted_tissue$tissuemodel,
#'   bulkgroups = bulkgroups
#' )
#'
#' differential_densities <- plot_differential_densities(
#'   density,
#'   scembedding,
#'   groupA = "bgroup_1",
#'   groupB = "bgroup_2",
#'   embedding_x = "UMAP1",
#'   embedding_y = "UMAP2",
#'   colourby = "c.type"
#' )
#'
#' print(differential_densities)
#'
#' @export
#'
plot_differential_densities <- function(
  densitydata,
  scembedding = NULL,
  groupA = NA,
  groupB = NA,
  embedding_x = "UMAP1",
  embedding_y = "UMAP2",
  colourby = NULL
) {

  if (
    !tibble::is_tibble(densitydata) ||
      !all(c("x", "y", "z", "group") %in% names(densitydata))
  ) {
    stop("invalid densitydata. must have columns x, y, z and groups.")
  }
  # prepare differential density plot for two groups
  diffdensity <- plot_element_differential_density(
    densitydata,
    groupA = groupA,
    groupB = groupB,
    x = "x",
    y = "y",
    z = "z",
    group = "group"
  )

  # prepare embedding plot of UMAP embedded single cells
  # according to the single cell embedding
  if (!is.null(scembedding)) {
    embedding <- plot_element_embedding(
      scembedding, embedding_x, embedding_y, colourby
    )
  } else {
    embedding <- NULL
  }

  # produce one ggplot object containing all three plots
  return(
    ggplot2::ggplot() +
      diffdensity +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
      embedding
  )

}
