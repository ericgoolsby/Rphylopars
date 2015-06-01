\name{phylopars.predict}
\alias{phylopars.predict}
\title{
Ancestral state reconstruction and missing value tip reconstruction
}
\description{
This function performs ancestral state reconstruction (and missing value tip reconstruction) given a fitted \code{phylopars} object.
}
\usage{
phylopars.predict(PPE, nodes, tips, verbose = FALSE)
}
\arguments{
  \item{PPE}{
Fitted \code{phylopars} object
}
  \item{nodes}{
Optional argument: vector of which internal nodes to reconstruct.
}
  \item{tips}{
Optional argument: vector of which species values to predict.
}
  \item{verbose}{
Whether or not to print progress of node reconstruction / tip prediction (default is FALSE).
}
}
\value{
Returns a list of reconstructed node and tip means, as well as standard deviations for each mean.
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

phylopars.predict(PPE)
}