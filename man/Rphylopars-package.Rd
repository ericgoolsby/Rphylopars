\name{Rphylopars-package}
\alias{Rphylopars-package}
\docType{package}
\title{
Phylogenetic Comparative Tools for Missing Data and Within-Species Variation
}
\description{
Tools for performing phylogenetic comparative methods for datasets with with multiple observations per species (intraspecific variation or measurement error) and/or missing data. Performs ancestral state reconstruction and missing data imputation on the estimated evolutionary model, which can be specified as Brownian Motion, Ornstein-Uhlenbeck, Early-Burst, Pagel's lambda, kappa, or delta, or a star phylogeny.
}
\details{
\tabular{ll}{
Package: \tab Rphylopars\cr
Type: \tab Package\cr
Version: \tab 0.2.3\cr
Date: \tab 2016-03-25\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Eric W. Goolsby, Jorn Bruggeman, Cecile Ane

Maintainer: Eric W. Goolsby \email{eric.goolsby.evolution@gmail.com}
}
\references{
Bruggeman J, Heringa J and Brandt BW. (2009) PhyloPars: estimation of missing parameter values using phylogeny. \emph{Nucleic Acids Research} 37: W179-W184.

Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
}
\keyword{ Phylogenetic variance-covariance }
\keyword{ Ancestral state reconstruction }
\keyword{ Intraspecific variation }
\keyword{ Within-species variation }
\keyword{ Measurement error }
\keyword{ Missing data }
\keyword{ Phylogenetic imputation }
\keyword{ package }
\examples{
# simulate data
sim_data <- simtraits(ntaxa = 15,ntraits = 4,nreps = 3,nmissing = 10)

# estimate parameters under Brownian motion
# pheno_error = TRUE assumes intraspecific variation
# pheno_correlated = FALSE assumes intraspecific variation is not correlated
# phylo_correlated = TRUE assumed traits are correlated

PPE <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree,
  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)

PPE

PPE$anc_recon # Ancestral state reconstruction and species mean prediction
PPE$anc_var # Prediction variance


###NOT RUN
# estimate parameters under multivariate OU
# PPE_OU <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree,
#    model="mvOU",pheno_error = TRUE,phylo_correlated = TRUE,
#    pheno_correlated = TRUE)
#
# PPE

}