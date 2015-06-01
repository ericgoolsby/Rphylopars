\name{phylopars.pca}
\alias{phylopars.pca}
\title{
Phylogenetic principal components analysis
}
\description{
This function, adapted from \code{\link[phytools]{phyl.pca}} in the \pkg{phytools} package, performs phylogenetic principal components analysis on a fitted \code{phylopars} object. Additionally, fast algorithms are used to avoid large matrix inversions using the \code{\link[phylolm]{three.point.compute}} function in \pkg{phylolm}. 
}
\usage{
phylopars.pca(PPE, mode = "cov", pgls_means = FALSE)
}
\arguments{
  \item{PPE}{
Fitted phylopars object
}
  \item{mode}{
Either "cov" (default) or "corr"
}
  \item{pgls_means}{
Whether to use raw species trait means of PGLS species trait means in forming the PCA
}
}
\value{
Returns an object of class \code{phyl.pca}.
}
\references{
Bruggeman J, Heringa J and Brandt BW. (2009) PhyloPars: estimation of missing parameter values using phylogeny. \emph{Nucleic Acids Research} 37: W179-W184.

Revell, L. J. (2009) Size-correction and principal components for interspecific comparative studies. \emph{Evolution}, 63, 3258-3268.

Ho, L. S. T. and Ane, C. (2014) A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. \emph{Systematic Biology}, 63: 397-408.
}
\author{
Eric W. Goolsby \email{eric.goolsby.evolution@gmail.com}, Cecile Ane, Jorn Bruggeman
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

# perform phylogenetic PCA
PPE_pca <- phylopars.pca(PPE)
summary(PPE_pca)
biplot(PPE_pca)
}