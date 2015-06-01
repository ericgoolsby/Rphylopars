\name{phylopars.rerun}
\alias{phylopars.rerun}
\title{
Rerun phylopars optimization
}
\description{
This function reruns parameter opimization for a phylopars object if convergence issues are suspected.
}
\usage{
phylopars.rerun(PPE)
}
\arguments{
  \item{PPE}{
Fitted \code{phylopars} object
}
}
\value{
Returns a fitted object of class \code{phylopars}.
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

# suspect convergence failure, so optimize parameters again
PPE <- phylopars.rerun(PPE)
}