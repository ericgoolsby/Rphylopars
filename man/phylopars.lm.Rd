\name{phylopars.lm}
\alias{phylopars.lm}
\title{
Phylopars regression
}
\description{
Performs phylogenetic regression.
}
\usage{
phylopars.lm(formula, trait_data, tree, model = "BM", pheno_error = TRUE,
phylo_correlated = TRUE, pheno_correlated = FALSE, calc_pheno = FALSE,
calc_pheno_auto_n = 20, calc_pheno_phenocovs, use_means = FALSE,
species_identifier = "species", verbose = FALSE, phylocov_start,
phenocov_start, theta_start, model_start, skip_optim = FALSE,
REML = TRUE, optim_limit = 50, BM_first = TRUE, usezscores = TRUE)
}
\arguments{
  \item{formula}{
A model formula.
}
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
A fitted phylopars.lm object.
}
\examples{

# no missing data
s <- simtraits(ntaxa = 50,ntraits = 3,nreps = 1,nmissing = 0)
p_BM <- phylopars.lm(V3~V1+V2,trait_data = s$trait_data,tree = s$tree,model = "BM")
p_lambda <- phylopars.lm(V3~V1+V2,trait_data = s$trait_data,tree = s$tree,model = "lambda")
summary(p_BM)
summary(p_lambda)

# with missing data
# notice that logLik, AIC, and BIC are not reported if data is missing
s <- simtraits(ntaxa = 50,ntraits = 3,nreps = 1,nmissing = 15)
p_BM <- phylopars.lm(V3~V1+V2,trait_data = s$trait_data,tree = s$tree,model = "BM")
p_lambda <- phylopars.lm(V3~V1+V2,trait_data = s$trait_data,tree = s$tree,model = "lambda")
summary(p_BM)
summary(p_lambda)
}