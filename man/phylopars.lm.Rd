\name{phylopars.lm}
\alias{phylopars.lm}
\title{
Rphylopars regression
}
\description{
Performs phylogenetic regression.
}
\usage{
phylopars.lm(formula, trait_data, tree, model = "BM", pheno_error,
phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE,
full_alpha = TRUE, phylocov_start, phenocov_start, model_par_start,
phylocov_fixed, phenocov_fixed, model_par_fixed, skip_optim = FALSE,
skip_EM = FALSE, EM_Fels_limit = 1000, repeat_optim_limit = 1,
EM_missing_limit = 50,repeat_optim_tol = 0.01, model_par_evals = 10,
max_delta = 10000, EM_verbose = FALSE,optim_verbose = FALSE, npd = FALSE,
nested_optim = FALSE, usezscores = TRUE, phenocov_list = list(),
ret_args = FALSE, ret_level = 1)
}
\arguments{
  \item{formula}{
Model formula -- e.g. Y~X1+X2
}
  \item{trait_data}{
A data frame with the first column labeled "species" (with species names matching tips on the phylogeny) and one column per trait. Each row corresponds to a single observation, and multiple observations for species are allowed. Missing data should be represented with NA.
}
  \item{tree}{
An object of class \code{phylo}
}
  \item{model}{
Model of evolution. Default is "BM". Alternative evolutionary models include "mvOU" (for the multivariate Ornstein-Uhlenbeck), or univariate tree transformations: "OU" "lambda", "kappa", "delta", "EB", "star".
}
  \item{pheno_error}{
If TRUE (default, unless <=1 observation per species is provided), parameters are estimated assuming within-species variation.
}
  \item{phylo_correlated}{
If TRUE (default), parameters are estimated assuming traits are correlated.
}
  \item{pheno_correlated}{
If TRUE (default), parameters are estimated assuming within-species observations traits are correlated.
}
  \item{REML}{
If TRUE (default), the algorithm will return REML estimates. If FALSE, maximum likelihood estimates will be returned.
}
  \item{full_alpha}{
Only applicable for the multivariate OU (model="mvOU"). If TRUE (default), a fully parametrized alpha matrix is fit. If FALSE, a diagonal alpha matrix is fit.
}
  \item{phylocov_start}{
Optional starting value for phylogenetic trait variance-covariance matrix. Must be of dimension n_traits by n_traits.
}
  \item{phenocov_start}{
Optional starting value for phenotypic trait variance-covariance matrix. Must be of dimension n_traits by n_traits.
}
  \item{model_par_start}{
Optional starting parameters for the evolutionary model. For model="mvOU", must be of dimension n_traits by n_traits. Otherwise, must be a single value.
}
  \item{phylocov_fixed}{
Optional fixed value for phylogenetic trait variance-covariance matrix. Must be of dimension n_traits by n_traits.
}
  \item{phenocov_fixed}{
Optional starting value for phenotypic trait variance-covariance matrix. Must be of dimension n_traits by n_traits.
}
  \item{model_par_fixed}{
Optional fixed parameter for the evolutionary model. For model="mvOU", must be of dimension n_traits by n_traits. Otherwise, must be a single value.
}
  \item{skip_optim}{
Whether to skip BFGS optimization (not recommended unless all parameters are fixed).
}
  \item{skip_EM}{
Whether to skip Expectation-Maximiation prior to generating starting parameters for BFGS optimization (not recommended unless providing fixed parameters).
}
  \item{EM_Fels_limit}{
Whether to skip Expectation-Maximiation prior to generating starting parameters for BFGS optimization (not recommended unless providing fixed parameters).
}
  \item{repeat_optim_limit}{
The number of times to repeat numerical optimization (default is 1).
}
  \item{EM_missing_limit}{
Maximum number of iterations for EM.
}
  \item{repeat_optim_tol}{
Maximum tolerance for repeated numerical optimization (only relevant if repeat_optim_limit>1).
}
  \item{model_par_evals}{
Number of times to evaluate univariate tree transformation models along the range of possible parameter values. Used to generate informed starting values for alternative evolutionary models if nested_optim=TRUE.
}
  \item{max_delta}{
Maximum allowed difference between the log-likelihood for EM-generated starting parameters and new parameters tried under numerical optimization. Extremely large deltas are likely to be numerical artifiacts. Prevents artificial convergence.
}
  \item{EM_verbose}{
Whether to print the log-likelihood during Expectation-Maximization.
}
  \item{optim_verbose}{
Whether to print log-likelihooods during numerical optimization.
}
  \item{npd}{
Whether to find the nearest positive-definite matrix for all covariance matrices during numerical optimization (slow -- only set to TRUE if converging to singular matrices).
}
  \item{nested_optim}{
Only relevant if fitting a univariate alternative evolutionary model. Tries multiple tree transformation parameter values along the range of possible values to make informed starting parameters. Slower than the default (nested_optim=FALSE), in which all parameters are estimated simultaneously.
}
  \item{usezscores}{
Whether or not ot use centered and standardized data during numerical optimization (recommended).
}
  \item{phenocov_list}{
An optional named list of species-specific within-species covariance matrices to be held fixed, as in Ives et al (2007). This option forces pheno_error and pheno_correlated to be FALSE, and uses mean species values instead of raw data. Raw variance should be divided by the number of observations per species (i.e., squared standard errors). See Ives et al (2007) for more details.
}
  \item{ret_args}{
For internal use only.
  \item{ret_level}{
For internal use only.
}
}
\value{
A fitted phylopars.lm object.
}
\examples{
# simulate data
sim_data <- simtraits(ntaxa = 15,ntraits = 4)

phylopars.lm(V4~V1+V2+V3,trait_data=sim_data$trait_data,tree=sim_data$tree)
}
