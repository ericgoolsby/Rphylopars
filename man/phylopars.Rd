\name{phylopars}
\alias{phylopars}
\title{
Estimation of phylogenetic and phenotypic covariance parameters
}
\description{
This function estimates parameters for the phylogenetic and phenotypic variance-covariance matrices for datasets with missing observations and multiple within-species observations. This function can also be used to fit altnerative evolutionary models, including Ornstein-Uhlenbeck, Early-Burst, star phylogeny, or Pagel's lambda, kappa, or delta. Reconstructed ancestral states and predicted species means (i.e., for missing data), along with prediction variances, are also provided.
}
\usage{
phylopars(trait_data, tree, model = "BM", pheno_error, phylo_correlated = TRUE,
pheno_correlated = TRUE, REML = TRUE, full_alpha = TRUE, phylocov_start,
phenocov_start, model_par_start, phylocov_fixed, phenocov_fixed, model_par_fixed,
skip_optim = FALSE, skip_EM = FALSE, EM_Fels_limit = 1000, repeat_optim_limit = 1,
EM_missing_limit = 50, repeat_optim_tol = 0.01, model_par_evals = 10, max_delta = 10000,
EM_verbose = FALSE, optim_verbose = FALSE, npd = FALSE,
nested_optim = FALSE, usezscores = TRUE)
}
\arguments{
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
}
\value{
An object of class \code{phylopars}. For models with phenotypic (within-species) covariance, the estimated percentage of variance explained by the phylogeny is provided as 100*(1 - phenotypic_variance/raw_variance), where raw_variance is the variance of all observations for a given trait across species (\code{var(PPE$trait_data[,2:ncol(PPE$trait_data)],na.rm=TRUE)}).
\item{logLik }{The log-likelihood of the model}
\item{pars }{A list composed of phylogenetic trait covariance and phenotypic (within-species) trait covariance, if estimated}
\item{model }{The model of evolution (e.g., BM, OU, lambda, etc.), and any additional evolutionary model parameters estimated. For OU models, stationary covariance is calculated from both phylogenetic covariance (Sigma) and alpha (see Supplement 1 of Clavel et al. 2015).}
\item{mu }{The estimate ancestral state at the root of the tree.}
\item{npars }{The total number of parameters estimated by optimization (used for AIC and BIC).}
\item{anc_recon }{Reconstructed ancestral states and species means. Row names correspond to species names (for the first 1:nspecies rows), and the remaining row names correspond to node numbers on a tree with edges in postorder: \code{reorder(tree,"postorder"). Or, if node labels were included on the original tree, row names correspond to node labels.}}
\item{anc_var }{Variance of reconstructed ancestral estimates and imputed species means.}
\item{anc_cov }{Covariance of estimates among variables.}
\item{tree }{The phylogenetic tree supplied to \code{phylopars}}
\item{trait_data }{The trait data supplied to \code{phylopars}}
\item{REML }{\code{TRUE} if REML, FALSE if ML.}
}
\references{
Bruggeman J, Heringa J and Brandt BW. (2009) PhyloPars: estimation of missing parameter values using phylogeny. \emph{Nucleic Acids Research} 37: W179-W184.

Clavel, J., Escarguel, G. & Merceron, G. (2015) mvmorph: an r package for fitting multivariate
261 evolutionary models to morphometric data. Methods in Ecology and Evolution, 6, 131-1319.

Felsenstein, J. (2008) Comparative methods with sampling error and within-species variation: contrasts revisited and revised. American Naturalist, 171, 713-725.

Ho L.S.T., Ane C. 2014. A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. Syst. Biol. 63:397-408.
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
