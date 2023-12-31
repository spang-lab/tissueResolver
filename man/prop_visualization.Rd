% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop_visualization.R
\name{prop_visualization}
\alias{prop_visualization}
\title{plots relative cell type fraction ratios between two groups}
\usage{
prop_visualization(props, grouping, groupA, groupB)
}
\arguments{
\item{props}{dataframe containing cell type proportions (output of cell_proportions)}

\item{grouping}{a dataframe containing "bulk_id" and "group" columns for group association of samples}

\item{groupA}{group that is in the numerator of the ratio}

\item{groupB}{group that is in the denominator of the ratio}
}
\value{
a ggplot object
}
\description{
plots relative cell type fraction ratios between two groups
}
