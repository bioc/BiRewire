\encoding{UTF-8}
\name{birewire.bipartite.from.incidence}

\alias{birewire.bipartite.from.incidence}

\title{Converts an incidence matrix into a bipartite graph. }
\description{This function creates an \code{\link{igraph}} bipartite graph from an incidence matrix.}
\usage{birewire.bipartite.from.incidence(matrix,directed=FALSE)}
\arguments{
  \item{matrix}{incidence matrix: an (n-by-m) binary matrix where rows correspond to vertices in the frist class while 
  columns correspond to vertices in the second one;}
  \item{directed}{Logical, if TRUE a directed graph is created.}
} 
\details{
The function calls \code{\link{graph.incidence}} of package \code{\link{igraph}}. See \code{\link{igraph}} documentation for more details.

}
\value{
Bipartite \emph{igraph} graph.

}
\references{
 Csardi, G. and  Nepusz, T (2006)
\emph{The igraph software package for complex network research},
 InterJournal, Complex Systems 
    url {http://igraph.sf.net}
  }
\author{Andrea Gobbi\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>
}
\examples{

library(igraph)
library(BiRewire)
g <-  graph.bipartite( rep(0:1,length=10), c(1:10))

##gets the incidence matrix of g
 m<-as.matrix(get.incidence(graph=g))

##rewire the current graph 
m2=birewire.rewire.bipartite(m,100)

#create the rewired bipartite graph
g2<-birewire.bipartite.from.incidence(m2,TRUE)
}
\keyword{bipartite graph, incidece matrix}
