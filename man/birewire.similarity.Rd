\encoding{UTF-8}
\name{birewire.similarity}

\alias{birewire.similarity}

\title{Compute the Jaccard similarity index between two binary matrices with the same number of non-null entries and the sam row- and column-wise sums.}
\description{Compute the Jaccard similarity index between two binary matrices with the same number of non-null entries and the sam row- and column-wise sums.} 
\usage{birewire.similarity( m1,m2)}
\arguments{
  \item{m1}{First matrix;}
  \item{m2}{Second matrix.  }
  }


\details{
The \strong{Jaccard} index between two sets \emph{M} and \emph{N} is defined as:

\eqn{{|M \cup N|}/{|M \cap N  |}   } 

 With M and N binary matrices, the Jaccard index is computed as:
\deqn{\frac{\sum N_{i,j} \wedge M_{i,j}}{ \sum N_{i,j} \vee M_{i,j}}.}\cr


The Jaccard index ranges between 0 and 1.


}

\value{
Returns the Jaccard similarity index between the two matrices
}

\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>
}
\examples{

library(igraph)
library(BiRewire)
g <- simplify(graph.bipartite( rep(0:1,length=100),
c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
seq(1,100,17),seq(1,100,19),seq(1,100,23),100
)))
g2=birewire.rewire.bipartite(g) 

birewire.similarity(get.incidence(g,sparse=FALSE),get.incidence(g2,sparse=FALSE))

}

\keyword{similarity,jaccard}
