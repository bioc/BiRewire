\encoding{UTF-8}
\name{birewire.bipartite.from.incidence}

\alias{birewire.bipartite.from.incidence}

\title{Converts an incidence matrix into a bipartite graph. }
\description{This function creates an \code{igraph} bipartite graph from an incidence matrix.}
\usage{birewire.bipartite.from.incidence(matrix,directed=FALSE,reverse=FALSE)}
\arguments{
  \item{matrix}{incidence matrix: an (n-by-m) binary matrix where rows correspond to vertices in the frist class while 
  columns correspond to vertices in the second one;}
  \item{directed}{Logical, if TRUE a directed graph is created.}
  \item{reverse}{Logical, if TRUE the edges will be directed from the nodes in the second class to those in the first one.}
} 
\details{
The function calls \code{graph.bipartite} of package \code{igraph}. See \code{igraph} documentation for more details.

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
g <- simplify(graph.bipartite( rep(0:1,length=100),
c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
seq(1,100,17),seq(1,100,19),seq(1,100,23),100
)))

##gets the incidence matrix of g
 m<-as.matrix(get.incidence(graph=g))

##rewire the current graph 
m2=birewire.rewire.bipartite(m,100)

#create the rewired bipartite graph
g2<-birewire.bipartite.from.incidence(m2,TRUE,FALSE)
}
\keyword{bipartite graph, incidece matrix}
