# This script generates reference output files for tissueResolver package tests.
# It runs key functions from the tissueResolver package on test data and saves
# the outputs using the version 1.0.0 of the package.
# If the package is updated, this script can be rerun to generate new reference
# files. Just uncomment all lines and run the script.

# ################################################################################
# ############################# load packages ####################################
# ################################################################################

# # load original tissueResolver package
# library(tissueResolver) # This should be version 1.0.0
# # load other required packages
# library(dplyr)
# library(purrr)
# library(ggplot2)

# # output path
# out_path <- "tests/testdata" # wd should be tissueResolver root directory

# ############



# ################################################################################
# ############################### fit_tissue #####################################
# ################################################################################

# # test input data
# test_bulkdata <- matrix(
#   c(
#     5, 3, 0, 1,
#     4, 0, 2, 1,
#     0, 2, 3, 5,
#     0, 1, 4, 4,
#     1, 0, 5, 3
#   ),
#   nrow = 5,
#   ncol = 4,
#   dimnames = list(
#     paste0("Gene", 1:5),
#     paste0("Bulk", 1:4)
#   )
# )

# test_sclibrary <- matrix(
#   c(
#     1, 0, 0, 0, 0,
#     1, 1, 0, 0, 0,
#     0, 1, 1, 0, 0,
#     0, 0, 1, 1, 0,
#     0, 0, 0, 1, 1
#   ),
#   nrow = 5,
#   ncol = 5,
#   dimnames = list(
#     paste0("Gene", 1:5),
#     paste0("Cell", 1:5)
#   )
# )


# # call fit_tissue_noboot function from tissueResolver package
# out_fit_tissue_noboot <- tissueResolver:::fit_tissue_noboot(
#   test_bulkdata,
#   test_sclibrary
# )

# # add bootstrap = FALSE to output (added in new function)
# out_fit_tissue_noboot$bootstrap <- FALSE

# # add convergence scores
# out_fit_tissue_noboot$log_store <- out_fit_tissue_noboot$log_store %>%
#   dplyr::mutate(
#     convergence = 0
#   )

# saveRDS(
#   out_fit_tissue_noboot,
#   file = file.path(out_path, "out_fit_tissue_noboot.rds")
# )

# ############



# ################################################################################
# ########################### specific_expression ################################
# ################################################################################

# # define inputs
# mapping <- dplyr::tibble(
#   cell_id = colnames(test_sclibrary),
#   celltype = c("CTA", "CTB", "CTA", "CTB", "CTA")
# )

# # fast_specific_expression_regulation
# out_fast_csre_noboot <- tissueResolver:::fast_specific_expression_regulation(
#   out_fit_tissue_noboot,
#   test_sclibrary,
#   mapping
# )

# # specific_expression_regulation
# out_csre_noboot <- tissueResolver:::specific_expression_regulation(
#   out_fit_tissue_noboot,
#   test_sclibrary,
#   mapping
# )

# out_csre <- tissueResolver:::specific_expression_regulation(
#   out_fit_tissue,
#   test_sclibrary,
#   mapping
# )

# # function to reorder output dataframes
# fix_csre <- function(csre) {

#   if ("regulation_boot" %in% colnames(csre)) {
#     csre <- csre %>%
#       dplyr::mutate(
#         regulation_boot = list(
#           purrr::map(
#             regulation_boot,
#             ~ifelse(is.nan(.x), 0, .x)) %>% unlist()
#         ),
#         regulation = ifelse(
#           is.nan(regulation),
#           mean(unlist(regulation_boot)),
#           regulation
#         )
#       )
#   }

#   # reorder csre to match new function output
#   csre <- csre %>%
#     # oder by bulk_id, celltype, gene
#     dplyr::arrange(bulk_id, celltype, gene) %>%
#     # move total_explained to the end
#     dplyr::filter(celltype != "total_explained") %>%
#     dplyr::bind_rows(
#       csre %>%
#         dplyr::arrange(bulk_id, celltype, gene) %>%
#         dplyr::filter(celltype == "total_explained")
#     ) %>%
#     dplyr::ungroup()

#   # remove rowwise_df class
#   class(csre) <- setdiff(class(csre), "rowwise_df")

#   csre
# }

# # reorder outputs to match
# out_fast_csre_noboot <- fix_csre(out_fast_csre_noboot)
# out_csre_noboot <- fix_csre(out_csre_noboot)
# out_csre <- fix_csre(out_csre)

# # save outputs
# saveRDS(
#   out_fast_csre_noboot,
#   file = file.path(out_path, "out_fast_csre_noboot.rds")
# )
# saveRDS(
#   out_csre_noboot,
#   file = file.path(out_path, "out_csre_noboot.rds")
# )
# saveRDS(
#   out_csre,
#   file = file.path(out_path, "out_csre.rds")
# )

# ############



# ################################################################################
# ############################ cell proportions ##################################
# ################################################################################

# # get cell proportions for original function
# out_prop_single <- tissueResolver:::cell_proportions_singlemodel(
#   out_fit_tissue_noboot,
#   mapping
# )

# out_prop_noboot <- tissueResolver::cell_proportions(
#   out_fit_tissue_noboot,
#   mapping
# )

# out_prop <- tissueResolver::cell_proportions(
#   out_fit_tissue,
#   mapping
# )

# # FUN to add percentage column (added in new function)
# add_per <- function(props, by = "celltype") {
#   # calculate cell type proportions/fractions
#   props <- props %>%
#     dplyr::group_by(bulk_id) %>%
#     dplyr::mutate(bulk_weight = sum(weight)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(per_cells = (weight / bulk_weight) * 100) %>%
#     dplyr::select(!bulk_weight) %>%
#     dplyr::rename(!!paste0(by, "_percent") := per_cells)
#   return(props)
# }

# # add percentage columns
# out_prop_noboot <- add_per(out_prop_noboot)
# out_prop <- add_per(out_prop)

# # save outputs
# saveRDS(
#   out_prop_single,
#   file = file.path(out_path, "out_prop_single.rds")
# )

# saveRDS(
#   out_prop_noboot,
#   file = file.path(out_path, "out_prop_noboot.rds")
# )

# saveRDS(
#   out_prop,
#   file = file.path(out_path, "out_prop.rds")
# )

# ############



# ################################################################################
# ############################## quality scores ##################################
# ################################################################################

# # fix function for no bootstrap case
# quality_scores_fix <- function(
#   csre,
#   bulks
# ) {

#   if (!"total_explained" %in% (csre %>% dplyr::pull(celltype) %>% unique())) {
#     stop(
#       "need total explained gene expression in celldf. 
#       Call specific_expression_regulation with compute_total=TRUE."
#     )
#   }

#   bulkdata <- tissueResolver:::promote_matrix(bulks)

#   if (ncol(bulkdata) < 20) {
#     warning("small number of bulks. quality scores will be unreliable.")
#   }

#   bulknames <- intersect(
#     colnames(bulkdata),
#     csre %>% dplyr::pull(bulk_id) %>% unique()
#   )

#   bulkdf <- tidyr::as_tibble(bulkdata) %>%
#     tibble::add_column(gene = rownames(bulkdata)) %>%
#     tidyr::pivot_longer(
#       -gene,
#       values_to = "bulk_expression",
#       names_to = "bulk_id"
#     )

#   total_fitted <- csre %>%
#     dplyr::filter(bulk_id %in% bulknames) %>%
#     dplyr::filter(celltype == "total_explained") %>%
#     dplyr::ungroup() %>%
#     dplyr::select(-celltype)

#   total_fitted <- total_fitted %>%
#     dplyr::inner_join(bulkdf, by=c("gene", "bulk_id"))

#   relres_score <- function(vec, bulkexpr) {
#     # this function works both for vectors as well as
#     # scalars and can be used in the bootstrap and the non-bootstrap version.
#     v <- abs(vec - bulkexpr) / (vec + bulkexpr)
#     # if denom is zero, the fit is still perfect w/o deviation
#     v[is.na(v)] <- 0.0
#     v
#   }

#   if ("expression_boot" %in% names(csre)) {

#     # is a bootstrap run.
#     res <- total_fitted %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(relres_boot = list(
#         relres_score(expression_boot, bulk_expression)
#       )) %>%
#       dplyr::mutate(relres = mean(relres_boot)) %>%
#       dplyr::mutate(relres_var = stats::var(relres_boot)) %>%  # relres
#       dplyr::ungroup() %>%
#       dplyr::select(
#         bulk_id, gene, relres, relres_var, expression, bulk_expression
#       )

#     gene_qc <- res %>%
#       dplyr::group_by(gene) %>%
#       dplyr::summarise(
#         relres_sample_var = stats::var(relres),
#         relres = mean(relres),
#         relres_mean_var = mean(relres_var),
#         expression = mean(expression),
#         bulk_expression = mean(bulk_expression)
#       )

#     bulk_qc <- res %>%
#       dplyr::group_by(bulk_id) %>%
#       dplyr::summarise(
#         relres_sample_var = stats::var(relres),
#         relres = mean(relres),
#         relres_mean_var = mean(relres_var)
#       )

#   } else {

#     # no bootstrap run
#     res <- total_fitted %>%
#       dplyr::mutate(relres =
#         relres_score(expression, bulk_expression)
#       ) %>%
#       dplyr::ungroup() %>%
#       dplyr::select(
#         bulk_id, gene, relres, expression, bulk_expression
#       )

#     gene_qc <- res %>% 
#       dplyr::group_by(gene) %>%
#       dplyr::summarise(
#         relres_sample_var = stats::var(relres),
#         relres = mean(relres),
#         expression = mean(expression),
#         bulk_expression = mean(bulk_expression)
#       )

#     bulk_qc <- res %>% 
#       dplyr::group_by(bulk_id) %>%
#       dplyr::summarise(
#         relres_sample_var = stats::var(relres),
#         relres = mean(relres)
#       )
#   }

#   list(genes = gene_qc, bulks = bulk_qc)

# }

# # get quality scores for original function
# out_qc_noboot <- quality_scores_fix(
#   out_csre_noboot,
#   test_bulkdata
# )

# out_qc <- tissueResolver::quality_scores(
#   out_csre,
#   test_bulkdata
# )

# # save outputs
# saveRDS(
#   out_qc_noboot,
#   file = file.path(out_path, "out_qc_noboot.rds")
# )

# saveRDS(
#   out_qc,
#   file = file.path(out_path, "out_qc.rds")
# )

# ############



# ################################################################################
# ########################### weight_visualization ###############################
# ################################################################################

# ### kde2d_weighted test ###
# # input data
# tmp_frame <- data.frame(
#   x = c(
#     rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)
#   ),
#   y = c(rep(x = 1:3, each = 3)),
#   c.type = c("A", "B", "C", "A", "C", "A", "B", "B", "C"),
#   group_0 = c(1, 0.5, 0.25,  rep(0, 6)),
#   group_1 = c(rep(0, 3), rep(1, 3), rep(0, 3)),
#   group_2 = rev(c(rep(1, 3), rep(0, 6)))
# )

# out_kde2d <- tissueResolver:::kde2d_weighted(
#   x = tmp_frame$x,
#   y = tmp_frame$y,
#   weight = tmp_frame$group_0
# )

# # save output
# saveRDS(
#   out_kde2d,
#   file = file.path(out_path, "out_kde2d.rds")
# )

# ###



# ### bulk_group_input test ###

# # Case 1
# vec_bulkgrp <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_3")

# bulk_group_input_1 <- tissueResolver:::bulk_group_input(
#   vec_bulkgrp,
#   bulkidcol = "bulk_id",
#   groupcol = "group"
# )

# # Case 2
# tbl_bulkgrp <- tibble(
#   "bulk_id" = c("b1", "b2", "b3"),
#   "group" = c("group_1", "group_2", "group_3")
# )

# bulk_group_input_2 <- tissueResolver:::bulk_group_input(
#   tbl_bulkgrp,
#   bulkidcol = "bulk_id",
#   groupcol = "group"
# )

# ###



# ### density_distribution test ###

# # define inputs
# sc_emb <- data.frame(
#   "cell_id" = paste0("cell_", 1:9),
#   "UMAP1" = rep(c(1, 2, 3), each = 3),
#   "UMAP2" = rep(x = 1:3, each = 3),
#   "celltype" =  c("A", "B", "C", "A", "C", "A", "B", "B", "C")
# )

# tm <- dplyr::tibble(
#   "fitted_genes" = NA,
#   "tissuemodel" = tibble(
#     "bulk_id" = c("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
#     "cell_id" = c(
#       paste0("cell_", 1:3), paste0("cell_", 4:6), paste0("cell_", 1:3)
#     ),
#     "weight" = c(1, 0.5, 0.25, 0.25, 0.5, 1, 0.3, 0.3, 0.3)
#   )
# )

# bulkgrp <- c("b1" = "group_1", "b2" = "group_2", "b3" = "group_1")

# # call function
# out_den_dist <- density_distribution(
#   scembedding = sc_emb,
#   tissuemodel = tm$tissuemodel,
#   bulkgroups = bulkgrp
# )

# # save output
# saveRDS(
#   out_den_dist,
#   file = file.path(out_path, "out_den_dist.rds")
# )

# ###



# ### plot_element_embedding test ###

# sc_emb <- data.frame(
#   "cell_id" = paste0("cell_", 1:9),
#   "UMAP1" = rep(c(1, 2, 3), each = 3),
#   "UMAP2" = rep(x = 1:3, each = 3),
#   "celltype" =  c("A", "B", "C", "A", "C", "A", "B", "B", "C")
# )

# embedding_plot <- tissueResolver:::plot_element_embedding(
#   scembedding = sc_emb,
#   colourby = "celltype"
# )

# ############



# ################################################################################
# ################################# utils ########################################
# ################################################################################

# # test promote_matrix function
# test_vec <- 1:3

# out_m <- tissueResolver:::promote_matrix(test_vec)

# test_m <- matrix(1:3, nrow = 3, ncol = 1)

# # test df_to_matrix function
# tbl <- tibble(
#   gene = paste0("Gene", 1:5) %>% rep(each = 3),
#   cell_id = paste0("Cell", 1:3) %>% rep(times = 5),
#   expression = c(
#     0, 1, 1, 0, 0,
#     1, 0, 1, 1, 0,
#     0, 0, 0, 1, 1
#   )
# )

# test <- tissueResolver:::df_to_matrix(tbl)

# test_rev <- tissueResolver:::matrix_to_df(test)

# all.equal(tbl, test_rev)

# ############