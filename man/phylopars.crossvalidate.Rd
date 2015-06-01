\name{phylopars.crossvalidate}
\alias{phylopars.crossvalidate}
\title{
Phylopars cross-validation
}
\description{
This function performs cross-validation of the fitted phylopars object using the leave-one-out method. For each species, all observations for a single trait are temporarily deleted, and the value is phylognetically imputed using phylopars.predict(). A regression between the observed and predicted species means is returned. Also estimated is mean bias and error of the phylopars model, as well as two alternative reference models: mean model and nearest neighbor model.
}
\usage{
phylopars.crossvalidate(PPE, plot = FALSE, verbose = FALSE)
}
\arguments{
  \item{PPE}{
Fitted \code{phylopars} object
}
  \item{plot}{
Whether or not to plot cross-validation statistics (default is TRUE)
}
  \item{verbose}{
Whether or not to print progress of cross-validation (default is FALSE)
}
}
\value{
Returns a summary of cross-validation statistics.
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

phylopars.crossvalidate(PPE)
}