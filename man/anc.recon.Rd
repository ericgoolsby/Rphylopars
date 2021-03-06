\name{anc.recon}
\alias{anc.recon}
\title{
Ultra-fast maximum likelihood ancestral state reconstruction
}
\description{
This function performs ancestral state reconstruction using a fast algorithm based on Ho and Ane (2014).
}
\usage{
anc.recon(trait_data, tree, vars = FALSE, CI = FALSE)
}
\arguments{
  \item{trait_data}{
A vector or matrix of trait values. Names or row names correspond to species names. Data cannot have any missing data or within-species variation (this type of data can be handled by the \link[Rphylopars]{phylopars} function).
}
  \item{tree}{
An object of class \code{phylo}.
}
  \item{vars}{
Whether to return the variances of the restricted maximum likelihood estimates
}
  \item{CI}{
Whether to return 95\% confidence intervals of the restricted maximum likelihood estimates
}
}
\value{
A named vector of maximum likelihood ancestral states (with names corresponding to node names if available or node numbers from the tree rearranged in postorder, as obtained by the command \code{reorder(tree,"postorder")}). If vars or CI is set to TRUE, a list is returned with these values included.
}
\author{
Felsenstein, J. (1985) Phylogenies and the comparative method. American Naturalist, 125, 1-15.

Ho L.S.T., Ane C. 2014. A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. Syst. Biol. 63:397-408.

Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.
}
\seealso{
\link[phytools]{fastAnc}, \link[ape]{ace}, \link[ape]{pic}
}
\examples{
require(ape)
tree <- rtree(10000) # random tree with 10,000 taxa
x <- setNames(rnorm(1e4),tree$tip.label) # random trait data
recon <- anc.recon(trait_data=x,tree=tree)
}