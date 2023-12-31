% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weight_visualization.R
\name{plot_element_embedding}
\alias{plot_element_embedding}
\title{plots the single cell embedding}
\usage{
plot_element_embedding(scembedding, x = "UMAP1", y = "UMAP2", colourby = NULL)
}
\arguments{
\item{scembedding}{A dataframe that holds 2D embeddings of the single cell data used to compute the tissue model}

\item{x}{column name of the x coordinate for the 2D embedding in scembedding}

\item{y}{column name of the y coordinate for the 2D embedding in scembedding}

\item{colourby}{if non-null, colour the single cells by this column of scembedding (e.g., celltype)}
}
\value{
a ggplot plot object that can further be customized, viewed or saved in a file. see ggplot2.
}
\description{
plots each UMAP point from the single cell embedding, optionally coloured by e.g., the associated cell type
}
\examples{
scembedding <- data.frame(
  "colnames" = paste0("cell_", 1:9),
  "UMAP1" = c(rep(x = 1, n = 3), rep(x = 2, n = 3), rep(x = 3, n = 3)),
  "UMAP2" = c(rep(x = 1:3, each = 3)),
  "c.type" = sample(c("CTA", "CTB", "CTC"), 9, replace = TRUE)
)

embedding_plot <- tissueResolver:::plot_element_embedding(
  scembedding,
  x = "UMAP1",
  y = "UMAP2",
  colourby = "c.type"
)

print(ggplot() +
  embedding_plot +
  coord_fixed())



}
