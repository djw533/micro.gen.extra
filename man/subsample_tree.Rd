% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggtree_function_extensions.R
\name{subsample_tree}
\alias{subsample_tree}
\title{Subsample a phylogenetic tree}
\usage{
subsample_tree(input_tree, cut_dist, plot_clades_to_prune = F)
}
\arguments{
\item{input_tree}{Tree object (e.g. from ape)}

\item{cut_dist}{max tip-root distance to select clades to be sub-sampled}

\item{plot_clades_to_prune}{plot original input tree with clades to be pruned annotated [Default = False]}
}
\value{
Pruned tree object
}
\description{
Subsample a phylogenetic tree by randomly selecting one node for all clades/subtrees of a maximum tip-root distance
}
\examples{
subsample_tree(input_tree = tree_object, cut_dist = 0.02)
}
