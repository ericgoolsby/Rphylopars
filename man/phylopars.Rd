\name{phylopars}
\alias{phylopars}
\title{
Estimation of phylogenetic and phenotypic covariance parameters
}
\description{
This function estimates parameters for the phylogenetic and phenotypic variance-covariance matrices for datasets with missing observations and multiple within-species observations. This function is a prerequisite for ancestral state reconstruction and missing tip estimation (\code{\link[Rphylopars]{phylopars.predict}}) and phylogenetic PCA (\code{\link[Rphylopars]{phylopars.pca}}).
}
\usage{
phylopars(trait_data, tree, model = "BM", pheno_error = TRUE,
phylo_correlated = TRUE, pheno_correlated = FALSE, calc_pheno = FALSE,
calc_pheno_auto_n = 20, calc_pheno_phenocovs, use_means=FALSE, 
species_identifier = "species", verbose = FALSE, phylocov_start,
phenocov_start, theta_start, model_start, skip_optim = FALSE, REML = TRUE,
optim_limit = 50, BM_first = TRUE, usezscores = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trait_data}{
A data frame with one column per trait, as well as a column labeled "species" (with species names matching tips on the phylogeny). Each row corresponds to a single observation, and multiple observation for species are allowed. Missing data should be represented with NA.
}
  \item{tree}{
An object of class \code{phylo}
}
  \item{model}{
Model of evolution. Default is "BM". Also accepts "lambda", "kappa", "delta", "EB", and "OUrandomRoot", "OUfixedRoot",
}
  \item{pheno_error}{
If TRUE (default), parameters are estimated assuming within-species variation.
}
  \item{phylo_correlated}{
If TRUE (default), parameters are estimated assuming traits are correlated.
}
  \item{pheno_correlated}{
If TRUE, parameters are estimated assuming within-species observations traits are correlated (default is FALSE).
}
  \item{calc_pheno}{
If TRUE, phenotypic error variance is estimated for each species. If the number of observations for a given species is less than calc_pheno_auto_n (default=20), pooled within-species variance is used, unless phenotypic covariance matrices are supplied in the calc_pheno_phenocovs argument.
}
  \item{calc_pheno_auto_n}{
The minimum number of within-species observations for a given species for which within-species variance is calculated individually for a given species rather than pooled within-species variance. Only relevant if calc_pheno=TRUE and calc_pheno_phenocovs is not supplied.
}
  \item{calc_pheno_phenocovs}{
User-supplied phenotypic covariance matrix OR a list of phenotypic covariance matrices for each species.
}
  \item{use_means}{
Whether to use means (as in Ives et al. 2007) or raw data. Default is raw data (use_means=FALSE).
}
  \item{species_identifier}{
Optional argument. The name of the column in \code{data} where species names are stored. The default value is "species", but can be changed for convenience if \code{data} has an alternative species column name.
}
  \item{verbose}{
Optional argument. Whether or not to print the log-likelihood during optimization (NOTE: constants in the likelihood calculation are dropped from the value reported during optimization.)
}
  \item{phylocov_start}{
Optional starting value for phylogenetic trait variance-covariance matrix. Must be of dimensions n_traits by n_traits.
}
  \item{phenocov_start}{
Optional starting value for phenotypic trait variance-covariance matrix. Must be of dimensions n_traits by n_traits.
}
  \item{theta_start}{
Optional starting parameters for \code{theta}.
}
  \item{model_start}{
Optional starting parameters for the evolutionary model.

}
  \item{skip_optim}{
Optional argument to skip optimization (NOTE: this will produce the WRONG answer unless BM=TRUE, there is no missing data, and pheno_error=FALSE; generally leave set to FALSE).
}
  \item{REML}{
If TRUE (default), the algorithm will return REML estimates. If FALSE, maximum likelihood estimates will be returned.

}
  \item{optim_limit}{
Optional: number of extra optimization attempts allowed if convergence is not achieved. Only necessary if convergence failure is suspected.
}
  \item{BM_first}{
Whether to estimate covariance assuming Brownian motion prior to estimating alternative evolutionary model paramters.
}
  \item{usezscores}{
Whether to convert data to z-scores during optimization (recommended).
}
}
\value{
An object of class \code{phylopars}.
}
\references{
Bruggeman J, Heringa J and Brandt BW. (2009) PhyloPars: estimation of missing parameter values using phylogeny. \emph{Nucleic Acids Research} 37: W179-W184.
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

# perform phylogenetic regression
p_BM <- phylopars.lm(V3~V1+V2,trait_data = sim_data$trait_data,tree = sim_data$tree,model = "BM")

# perform ancestral state reconstruction / tip prediction
PPP_ace <- phylopars.predict(PPE)
}
