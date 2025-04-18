% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xegaGpMutate.R
\name{xegaGpMutateFilterGene}
\alias{xegaGpMutateFilterGene}
\title{Mutate a gene (with a node filter)}
\usage{
xegaGpMutateFilterGene(g, lF)
}
\arguments{
\item{g}{Derivation tree.}

\item{lF}{Local configuration of the genetic algorithm.}
}
\value{
Derivation tree.
}
\description{
\code{xegaGpMutateGeneFilter()} replaces 
             a randomly selected subtree by
             a random derivation tree with the same root symbol 
             with a small probability.
             Only non-terminal nodes with a depth
             between \code{lF$MinMutInsertionDepth()} and
             \code{lF$MaxMutInsertionDepth()} are considered 
             as tree insertion points.
             Depth-bounds are respected.
}
\details{
Mutation is controlled by three local parameters: 
          \enumerate{
           \item \code{lF$MaxMutDepth()} controls the maximal depth of 
                 the new random generation tree.
           \item \code{lF$MinMutInsertionDepth()} and 
                 \code{lF$MaxMutInsertionDepth()} control the possible 
                 insertion points for the new random derivation tree.
                 The depth of the insertion node must be 
                 between \code{lF$MinMutInsertionDepth()} and
                 \code{lF$MaxMutInsertionDepth()}.
          }
}
\examples{
gene1<-xegaGpInitGene(lFxegaGpGene)
xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpMutateFilterGene(gene1, lFxegaGpGene)
xegaGpDecodeGene(gene, lFxegaGpGene)

}
\seealso{
Other Mutation: 
\code{\link{xegaGpMutateAllGene}()}
}
\concept{Mutation}
