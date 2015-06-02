\name{Rphylopars-package}
\alias{Rphylopars-package}
\docType{package}
\title{
Phylogenetic Comparative Tools for Missing Data and Within-Species Variation
}
\description{
Tools for performing phylogenetic comparative methods for datasets with with multiple observations per species (intraspecific variation or measurement error) and/or missing data. Available functions allow ancestral state reconstruction, missing data imputation, phylogenetic regression (phylogenetic generalized least squares), and phylogenetic principal components analysis (pPCA).
}
\details{
\tabular{ll}{
Package: \tab Rphylopars\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2015-06-02\cr
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
\keyword{ Phylogenetic regression }
\keyword{ Phylogenetic principle components analysis }
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
  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = FALSE)

PPE

# perform phylogenetic PCA
PPE_pca <- phylopars.pca(PPE)
summary(PPE_pca)

# perform phylogenetic regression
p_BM <- phylopars.lm(V3~V1+V2,trait_data = sim_data$trait_data,tree = sim_data$tree,model = "BM")

# perform ancestral state reconstruction / tip prediction
PPP_ace <- phylopars.predict(PPE)
}