\encoding{UTF-8}
\name{birewire.similarity.dsg}

\alias{birewire.similarity.dsg}
\title{Computes the Jaccard index between two directed signed networks.}
\description{Computes the Jaccard index between two directed signed networks, the two graphs must have the same dimensions.} 
\usage{birewire.similarity.dsg(m1,m2)}
\arguments{
 \item{m1}{A dsg object.}
  \item{m2}{A dsg object.}
 } 
\details{ Computes the Jaccard index between two directed signed networks, the two graphs must have the same dimensions.
}
\value{
A real number between 0 (totally dissimilar graphs) and 1 (equal graphs).}
\references{
 Iorio, F. and Gobbi, A.  and  Saez-Rodriguez, J. (2015) \emph{XXX} YYY.\cr
  }

\examples{

data(test_dsg)
dsg=birewire.induced.bipartite(test_dsg)
tmp= birewire.rewire.dsg(dsg,verbose=FALSE)
dsg2=birewire.build.dsg(tmp)
jaccard=birewire.similarity.dsg(dsg,tmp)
}
\keyword{directed graph, rewire, pathway, signaling}
