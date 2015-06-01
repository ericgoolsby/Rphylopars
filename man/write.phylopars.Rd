\name{write.phylopars}
\alias{write.phylopars}
\title{
Write data and tree files for Python phylopars compatability.
}
\description{
Writes data and tree files for Python phylopars compatibility.
}
\usage{
write.phylopars(trait_data, tree, data_file, tree_file, species_identifier = "species")
}
\arguments{
  \item{trait_data}{
A data frame with one column per trait, as well as a column labeled "species" (with species names matching tips on the phylogeny). Each row corresponds to a single observation, and multiple observation for species are allowed. Missing data should be represented with NA.
}
  \item{tree}{
An object of class \code{phylo}
}
  \item{data_file}{
Desired path to write data file.
}
  \item{tree_file}{
Desired path to write tree file.
}
  \item{species_identifier}{
Title of species column in data file. Defaulted to \code{"species"}
}
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

write.phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree,data_file = "data_path.txt",
tree_file = "tree_path.new")
}
