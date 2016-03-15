\encoding{UTF-8}
\name{birewire.similarity.dsg}

\alias{birewire.similarity.dsg}

\title{Compute the Jaccard similarity index between dsg.}
\description{Compute the Jaccard similarity index between dsg objects described in the same way (matrices of graphs).}  
\usage{birewire.similarity.dsg( m1,m2)}
\arguments{
  \item{m1}{First dsg;}
  \item{m2}{Second dsg.  }
  }


\details{
See \code{\link{birewire.similarity}} for more details. 

}

\value{
Returns the Jaccard similarity index between the objects.
}

\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>
}
\examples{

library(BiRewire)
data(test_dsg)
dsg <-  birewire.induced.bipartite(test_dsg,sparse=FALSE)
birewire.similarity.dsg(dsg,birewire.rewire.dsg(dsg))
dsg <-  birewire.induced.bipartite(test_dsg,sparse=TRUE)
birewire.similarity.dsg(dsg,birewire.rewire.dsg(dsg))
}

\keyword{similarity,jaccard}
