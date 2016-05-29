\name{fast.SSC}
\alias{fast.SSC}
\title{
Fast Phylogenetic Signal Using Sum of Squared Changes (SSC)
}
\description{
This function uses a fast ancestral state reconstruction algorithm (\link[Rphylopars]{anc.recon}, Goolsby, In review) to calculate the sum of squared changes bewteen ancestral and descendant nodes/tips, as described in Klingenberg and Gidaszewski (2010). Significance is assessed via phylogenetic permutation.
}
\usage{
fast.SSC(trait_data, tree, niter = 1000)
}
\arguments{
  \item{trait_data}{
A vector or matrix of trait values. Names or row names correspond to species names. Data cannot have any missing data or within-species variation.
}
  \item{tree}{
An object of class \code{phylo}.
}
  \item{niter}{
Number of iterations for hypothesis testing (default=1000).
}
}
\value{
\item{pvalue }{Description of 'comp1'}
\item{scaled.SSC }{Scaled sum of squared changes. A value less than 1 indicates less phylogenetic signal as measured by SSC than expected under Brownian motion, and a value greater than 1 indicates greater phylogenetic signal as measured by SSC than expected under Brownian motion.}
\item{SSC }{Total sum of squared changes (SSC)}
}
\references{
Goolsby E.W. (In Review). "Ultra-fast ancestral state reconstruction of continuous characters: a rerooting-free maximum likelihood approach."

Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. Systematic Biology. Accepted.

Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution, 57:717-745.

Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.

Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Systematic Biology. 63:685-697.
}
\author{
Eric W. Goolsby
}
\seealso{
\link[phytools]{phylosig}, \link[geomorph]{physignal}, \link[phylocurve]{fast.geomorph.physignal}
}
\examples{
sim_dat <- simtraits(ntaxa = 100,ntraits = 4)
fast.SSC(trait_data = sim_dat$trait_data,tree = sim_dat$tree)
}
