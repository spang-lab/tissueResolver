#' weighted 2D Kernel density estimation
#'
#' The following function is an adaption from MASS::kde2d, 
#'
#' @param x numeric, x coordinate of the embedding
#' @param y numeric, y coordinate of the embedding
#' @param weight numeric, weight for each data point
#' @param n.grid integer, number of grid points for plotting, defaults to 25
#' @param lims numeric with length of 4 c(lower(x),upper(x), lower(y), upper(y)), range of the plot, allowing for plotting only a subframe of the data.
#' Defaults to 'c(range(x), range(y))'
#' @param bandwidth bandwidth for the KDE, inducing the smoothing for the density estimate.
#' Numeric with length of 2,
#' defaults to 'c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))/4' 
#'
#' @return data.frame with 3 entries: 'x' and 'y' grid coordinates
#' in long format. 
#' 'z' (= weighted density) in long format
#' 
#' @examples
#' rm(list = ls()); gc()
#' library(tissueResolver)
#' library(ggplot2)
#' 
#' set.seed(1)
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
# 
kde2d_weighted <- function (
                      x,
                      y,
                      weight,
                      bandwidth = NULL,
                      n.grid = 25,
                      lims = c(range(x), range(y))
                  ) {
  # number of data points: 
  n.x <- length(x)
  
  if (length(y) != n.x) {
    stop("In 'kde2d_weighted': 'length(x)' and 'length(y)' differ.")
  }
  if (length(weight) != n.x) {
    stop("In 'kde2d_weighted': 'length(x)' and 'length(weight)' differ.")
  }
  
  grid.x <- seq(from = lims[1], to = lims[2], length = n.grid) 
  grid.y <- seq(from = lims[3], to = lims[4], length = n.grid) 
  
  if(is.null(bandwidth)){
    bandwidth <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))/4
  }
  
  # distance of each data point to each grid point scaled by the bandwidth,
  # for computation of kernel
  distance.x <- outer(grid.x, x, "-")/bandwidth[1] 
  distance.y <- outer(grid.y, y, "-")/bandwidth[2] 
  
  # weight matrix
  left.mat <- matrix(rep(weight, n.grid), nrow = n.grid, ncol = n.x, byrow = TRUE)

  # matrices for the kernel
  middle.mat <- matrix(dnorm(distance.x), n.grid, n.x)
  right.mat <- matrix(dnorm(distance.y), n.grid, n.x)
  
  # weighted kde differs from classical kde by multiplying the kernel with weight matrix
  # and dividing by the cumulated weights
  # for classical kde, see Section 5.6 in 
  # Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
  z <- (left.mat*middle.mat) %*% t(right.mat)/(sum(weight) * bandwidth[1] * bandwidth[2]) 
  
  # the returned dataframe stores the grid coordinates
  # and the weighted kde
  ret <- data.frame(
    x = rep(grid.x, n.grid)
    , y = rep(grid.y, each = n.grid)
    , z = as.vector(z)
  )
  return(ret)
}


#' input wrapper for accepting both named vectors as well as dataframes containing bulk groups
#' @param bulkgroups vector or dataframe holding the bulk-group association
#' @param bulkidcol if bulkgroups is a tibble, the name of the column that holds the bulk ids.
#' @param groupcol if bulkgroups is a tibble, the name of the column that holds the group ids.
#'
#' @return named vector of group ids. Names are the bulk ids, values are group ids.
#'
#' @examples
#' # Case 1
#' bulkgroups <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_3")
#' 
#' bulk_group_input_1 <- tissueResolver:::bulk_group_input(bulkgroups, bulkidcol = "bulk_id", groupcol = "group")
#' 
#' # Case 2
#' bulkgroups <- tibble(
#'   "bulk_id" = c("b1", "b2", "b3"),
#'   "group" = c("group_1", "group_2", "group_3")
#' )
#' 
#' bulk_group_input_2 <- tissueResolver:::bulk_group_input(bulkgroups, bulkidcol = "bulk_id", groupcol = "group")
#' 
#
bulk_group_input <- function(bulkgroups, bulkidcol = "bulk_id", groupcol="group") {
  if(tibble::is_tibble(bulkgroups)) {
    if( ! bulkidcol %in% names(bulkgroups)) {
      stop(paste0(bulkidcol, " is not a column of the bulk group definition dataframe"))
    }
    if( ! groupcol %in% names(bulkgroups)) {
      stop(paste0(groupcol, " is not a column of the bulk group definition dataframe"))
    }
    tmp <- bulkgroups[[groupcol]]
    names(tmp) <- bulkgroups[[bulkidcol]]
    return(tmp)
  } else if( is.vector(bulkgroups) ){
    if( is.null(names(bulkgroups)) ){
      stop("bulk group input vector is unnamed. Expect bulk ids as names of the group vector")
    }
    return(bulkgroups)
  } else {
    stop("invalid input for bulkgroups")
  }
}


#' density distributions of tissue models, allowing individual groups
#' @description computes weighted 2D-densities for cell populations fitted using fit_tissue.
#' If groups are given, such densities are computed for every group separately.
#'
#' @param scembedding A dataframe that holds 2D embeddings of the single cell data used to compute the tissue model
#' @param x name of the x coordinate of the 2D embedding.
#' @param y name of the y coordinate of the 2D embedding
#' @param scids name of the column holding the cell identifiers in the embedding data frame.
#' @param tissuemodel tissue model, output from fit_tissue
#' @param model_cellids name of the cell_id column in the model data frame
#' @param model_bulkids name of the bulk id column in the model data frame
#' @param bulkgroups an optional mapping of bulks to condition groups (e.g. treated/untreated, etc)
#' @param bulk_groupids names of the column that holds the group id within bulkgroups if that is a dataframe
#' @param verbose print more output
#'
#' @return dataframe with a kde density z-value for every x, y coordinate and an additional group column
#' 
#' @examples 
#' set.seed(1)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   # "colnames" = c(paste0("cell_", 1:3), paste0("cell_", 1:3), paste0("cell_", 1:3)),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#' 
#' fitted_tissue <- tibble(
#'   "fitted_genes" = NA,
#'   "tissuemodel" = tibble(
#'     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#'     # "cell_id" = paste0("cell_", 1:9),
#'     "cell_id" = c(paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)),
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
){

  # can deal with "full" tissuemodel, i.e. output of fit_tissue -> filter only tissuemodel tibble
  if( is_tissuemodel(tissuemodel) ) {
    tissuemodel <- tissuemodel$tissuemodel
  # can also deal with only tissuemodel tibble, then no filtering necessary 
  } else if( ! tibble::is_tibble(tissuemodel) ) {
    stop("tissuemodel is not a tissuemodel or dataframe")
  }
  
  if( ! model_bulkids %in% colnames(tissuemodel)) {
    stop(paste0("no bulk id group \"", model_bulkids, " in tissuemodel."))
  }
  if( is.null(bulkgroups) ) {
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
  if(! scids %in% colnames(scembedding)){
    stop(
      paste0(
        "In 'density_distribution': there is no '", scids, "' column, ",
        " in 'scembedding'. Therefore I can not map weights from ",
        "'tissuemodel' to scRNA-seq profiles"
      )
    )
  }
  
  if(!model_cellids %in% colnames(tissuemodel)){
    stop(
      paste0(
        "In 'density_distribution': there is no '", model_cellids, "' column, ",
        " in 'tissuemodel'. Therefore I can not map ",
        "weights from 'tissuemodel' to scRNA-seq profiles"
      )
    )
  }
  all.bulks <- names(bulkgroups)
  if(any(duplicated(all.bulks))){
    stop(
      paste0(
        "In 'density_distribution': 'names(bulkgroups)' must be unique."
      )
    )
  }
  # only consider unique bulk ids which both are given in the tissuemodel
  # and in the bulkgroups vector
  fitted.bulks <- unique(tissuemodel[[model_bulkids]])
  bulks.intersect <- intersect(all.bulks, fitted.bulks)

  if(length(bulks.intersect) < 2){
    stop(
      paste0(
        "In 'density_distribution': Only ",
        length(bulks.intersect), " bulks detected. Provide at least 2."
      )
    )
  }
  
  # end of sanity checks -------------------------------------------------------
  
  
  if(verbose){
    print(
      paste0(
        "--> Detected a bulk intersect of ", 
        length(bulks.intersect), ". Start to average per group."
      )
    )
  }
  # ATTENTION: Theoretically there could be different bulk ids in the
  # bulkgroup mapping compared to the tissuemodel's bulk ids
  unique.groups <- unique(bulkgroups)
  if(verbose){
    print(
      paste0(
        "--> Detected ", length(unique.groups), " groups." 
      )
    )
  }
  
  densities <- tibble()

  for(groupid in unique.groups){
    # filter the bulk ids from the boulkgroup mapping belonging to one group
    idx.this.groups <- which(bulkgroups == groupid)
    bulks.this.group <- names(bulkgroups)[idx.this.groups]
    
    # filter the bulk ids from the tissuemodel belonging to one bulk group
    idx.fitted.tissue <- which(
      tissuemodel[[model_bulkids]] %in% bulks.this.group
    )
    # store the part of the tissuemodel belonging to this group
    fits.for.this.group <- tissuemodel[idx.fitted.tissue, ]

    print(groupid)

    # stores for each cell in this group the mean across all weights
    # belonging to this cell
    weights.for.this.group <- sapply(
      X = unique(fits.for.this.group[[model_cellids]])
        # function that computes the mean over all weights
        # of a certain cell belonging to this group
      , FUN = function(cell){
        # retrieve the weight(s) for a certain cell within this group
        idx.for.this.cell <- which(fits.for.this.group[[model_cellids]] == cell)
        weights.for.this.cell <- as.numeric(
          fits.for.this.group[idx.for.this.cell, ][["weight"]]
        )
        return(
          mean(weights.for.this.cell)
        )
      }
    )
    weights.for.all.cells <- rep(0, nrow(scembedding))
    names(weights.for.all.cells) <- scembedding[[scids]]
    weights.for.all.cells[names(weights.for.this.group)] <- weights.for.this.group

    # store for each group the density matrix of the kde,
    # where the weight vector is given by the mean cell weight for each cell within this group
    densities <- rbind(densities,
      as_tibble(tissueResolver:::kde2d_weighted(
                         x = scembedding[[x]]
                       , y = scembedding[[y]]
                       , weight = weights.for.all.cells
                       )) %>% mutate(group = groupid))
  }
  return(densities)
}

#' plots the single cell embedding
#' @description plots each UMAP point from the single cell embedding, optionally coloured by e.g., the associated cell type
#' @param scembedding A dataframe that holds 2D embeddings of the single cell data used to compute the tissue model
#' @param x column name of the x coordinate for the 2D embedding in scembedding
#' @param y column name of the y coordinate for the 2D embedding in scembedding
#' @param colourby if non-null, colour the single cells by this column of scembedding (e.g., celltype)
#' 
#' @return a ggplot plot object that can further be customized, viewed or saved in a file. see ggplot2.
#'
#' @examples
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
#' 
#' 
plot_element_embedding <- function(scembedding, x="UMAP1", y="UMAP2",colourby=NULL) {
  if( is.null(colourby) ) {
    return(
      geom_point(data=scembedding, aes(x=.data[[x]], y=.data[[y]]), alpha=0.5)
    )
  } else {
    if(! colourby %in% names(scembedding)) {
      stop(paste0(colourby, " is not a column of scembedding."))
    }
    return(
      geom_point(data=scembedding, aes(x=.data[[x]], y=.data[[y]], colour=forcats::as_factor(.data[[colourby]])), alpha=0.5)
    )
  }
}

#' plots the differential density between two groups as a tile plot
#' @description having the kde density at each grid point of the single cell embedding for each bulkgroup provided by density_distribution, we
#' can assess the difference in densities between two of these groups. This differential densitiy is then visualized in the kde grid space as tile plot.
#' If group labels are omitted, they are tried to be guessed (if only two groups exist).
#' @param densities dataframe that is returned by density_distribution (see the help of that function)
#' @param groupA name of the group with the positive weights (appearing red in the plot where this group has a larger weight than groupB).
#' @param groupB name of the group with the negative relative weights (appearing blue in the plot where this group has a larger weight than groupA).
#' @param x column name of the x coordinate of the kde grid
#' @param y column name of the y coordinate of the kde grid
#' @param z column name of the z coordinate for the kde grid
#' @param group column name of the densities dataframe storing the bulkgroups
#' 
#' @return a ggplot plot object that can further be customized, viewed or saved in a file. see ggplot2.
#' @examples 
#' set.seed(1)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   # "colnames" = c(paste0("cell_", 1:3), paste0("cell_", 1:3), paste0("cell_", 1:3)),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#' fitted_tissue <- tibble(
#'   "fitted_genes" = NA,
#'   "tissuemodel" = tibble(
#'     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#'     # "cell_id" = paste0("cell_", 1:9),
#'     "cell_id" = c(paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)),
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
#' element_differential_density <- tissueResolver:::plot_element_differential_density(
#'   density,
#'   groupA = "bgroup_1",
#'   groupB = "bgroup_2",
#'   x = "x",
#'   y = "y",
#'   z = "z",
#'   group = "group"
#' )
#' 
#' print(ggplot() +
#'   element_differential_density +
#'   coord_fixed())

plot_element_differential_density <- function(
                                        densities,
                                        groupA = NA,
                                        groupB = NA,
                                        x="x",
                                        y="y",
                                        z="z",
                                        group="group") {
  groups <- unique(densities[[group]])
  ngroups <- length(groups)

  # if there are only two groups in the densities df and one of the input groups is NA,
  # specify the NA group as the remaining group
  if( ngroups == 2 ){
    if( is.na(groupA) && (! is.na(groupB)) ){
      groupA <- setdiff(groups, groupB)
    } else if(is.na(groupB) && (! is.na(groupA))) {
      groupB <- setdiff(groups, groupA)
    }
  } else if( is.na(groupA) || is.na(groupB) ){
    stop("no groups given. explicitly set groupA and groupB.")
  } else if( groupA == groupB ) {
    stop("two differential groups must not be equal.")
  }
  if( ( ! groupA %in% groups) || (! groupB %in% groups) ) {
    stop("groups for differential densities are not in the defined groups.")
  }
  # compute the differential densities for the two given groups
  diffdensity <- densities %>% 
    # filter only the rows of the df belonging to the two groups
    filter(.data[[group]] %in% c(groupA, groupB)) %>%
    # give each group a separate column
    pivot_wider(names_from=all_of(group), values_from=all_of(z), names_prefix="group_") %>%
    # give both groups the names specified by the input groups
    dplyr::rename(groupA = all_of(paste0("group_", groupA))) %>%
    dplyr::rename(groupB = all_of(paste0("group_", groupB))) %>%
    # store the differential density between groupA and groupB in the column diff
    mutate(diff = groupA - groupB) %>%
    # delete the seperate densities and keep only the differential density
    dplyr::select(all_of(c(x,y,"diff")))
  return(
    # represent the density as coloured tiles
    geom_tile(data=diffdensity, aes(x=.data[[x]], y=.data[[y]], fill=diff))
  )
}


#' plots differential group densities between two groups A and B.
#'
#' @description visualizes groups for which a tissue was fitted and their differences.
#' The difference in density between groupA and groupB is computed and plotted.
#' If group labels are omitted, they are tried to be guessed (if only two groups exist).
#' Besides the differential densities the plot contains the visualization of the single cell embedding
#' and a contour plot illustrating the density of each bulk group
#'
#' @param densitydata dataframe that is returned by density_distribution (see the help of that function)
#' @param scembedding embedding of the cells that the tissue model was computed on (dataframe or matrix-like). Optional. If NULL, no embedding will be plotted but just the densities.
#' @param groupA name of the group with the positive weights (appearing red in the plot where this group has a larger weight than groupB).
#' @param groupB name of the group with the negative relative weights (appearing blue in the plot where this group has a larger weight than groupA).
#' @param embedding.x column name of the x coordinate for the 2D embedding in scembedding
#' @param embedding.y column name of the y coordinate for the 2D embedding in scembedding
#' @param colourby if non-null, colour the single cells by this column of scembedding (e.g., celltype)
#'
#' @return a ggplot plot object that can further be customized, viewed or saved in a file. see ggplot2.
#' 
#' @examples
#' set.seed(1)
#' scembedding <- data.frame(
#'   "colnames" = paste0("cell_", 1:9),
#'   # "colnames" = c(paste0("cell_", 1:3), paste0("cell_", 1:3), paste0("cell_", 1:3)),
#'   "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
#'   "UMAP2" = c(rep(x = 1:3, each = 3)),
#'   "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
#' )
#' fitted_tissue <- tibble(
#'   "fitted_genes" = NA,
#'   "tissuemodel" = tibble(
#'     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#'     # "cell_id" = paste0("cell_", 1:9),
#'     "cell_id" = c(paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)),
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
#'   embedding.x = "UMAP1",
#'   embedding.y = "UMAP2",
#'   colourby = "c.type"
#' )
#' 
#' print(differential_densities)
#' @export
#'
#
plot_differential_densities <- function(
                                  densitydata,
                                  scembedding = NULL,
                                  groupA=NA,
                                  groupB=NA,
                                  embedding.x = "UMAP1",
                                  embedding.y = "UMAP2",
                                  colourby=NULL) {
  if( ! tibble::is_tibble(densitydata) || ! all(c("x", "y", "z", "group") %in% names(densitydata))) {
    stop("invalid densitydata. must have columns x, y, z and groups.")
  }
  # prepare differential density plot for two groups
  diffdensity <- plot_element_differential_density(
                                  densitydata,
                                  groupA = groupA,
                                  groupB = groupB,
                                  x="x",
                                  y="y",
                                  z="z",
                                  group="group")

  # prepare embedding plot of UMAP embedded single cells
  # according to the single cell embedding
  if( ! is.null(scembedding ) ) {
    embedding <- plot_element_embedding(scembedding, embedding.x, embedding.y, colourby)
  } else {
    embedding <- NULL
  }

  # produce one ggplot object containing all three plots
  return(
    ggplot() +
    diffdensity + scale_fill_gradient2(low="blue", mid="white", high="red") +
    embedding
  )
}
