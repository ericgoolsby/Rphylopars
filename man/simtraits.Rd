\name{simtraits}
\alias{simtraits}
\title{
Simulate traits for phylopars estimation
}
\description{
Simulates traits for code{phylopars} estimation.
}
\usage{
simtraits(ntaxa = 15, ntraits = 4, nreps = 3, nmissing = 10, tree, v, anc, intraspecific)
}
\arguments{
  \item{ntaxa}{
Either number of taxa (\code{ntaxa}) or a \code{tree} can be supplied.
}
  \item{ntraits}{
Number of traits to be simulated.
}
  \item{nreps}{
Number of replicates per trait per species to simulate.
}
  \item{nmissing}{
Number of randomly missing trait values.
}
  \item{tree}{
Either number of taxa (\code{ntaxa}) or a \code{tree} can be supplied.
}
  \item{v}{
Trait covariance (\code{v}) can be optionally supplied; otherwise off-diagonal elements are set to \code{0.8}.
}
  \item{anc}{
Value for ancestral state at root node for \code{\link[phytools]{fastBM}} from the \pkg{phytools} package.
}
  \item{intraspecific}{
Optional value for within-species variance.
}
}
\value{
\item{data}{Data for \code{phylopars()}}
\item{tree}{Tree for \code{phylopars()}}
}
\references{
Bruggeman J, Heringa J and Brandt BW. (2009) PhyloPars: estimation of missing parameter values using phylogeny. \emph{Nucleic Acids Research} 37: W179-W184.
}
\author{
Eric W. Goolsby \email{eric.goolsby.evolution@gmail.com}
}
\examples{
# simulate data
sim_data <- simtraits(ntaxa = 15,ntraits = 4,nreps = 3,nmissing = 10)

# estimate parameters under Brownian motion
# pheno_error = TRUE assumes intraspecific variation
# pheno_correlated = FALSE assumes intraspecific variation is not correlated
# phylo_correlated = TRUE assumed traits are correlated

PPE <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree,
  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = FALSE)

PPE
}