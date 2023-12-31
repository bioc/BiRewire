\encoding{UTF-8}
\name{birewire.rewire.undirected}

\alias{birewire.rewire.undirected}
\title{Efficient rewiring of undirected graphs}
\description{Optimal implementation of the switching algorithm. It returns the rewired version of the initial undirected graph or its adjacency matrix.} 
\usage{birewire.rewire.undirected(adjacency, max.iter="n",accuracy=0.00005,
verbose=TRUE,MAXITER_MUL=10,exact=FALSE)}
\arguments{
  \item{adjacency}{An \code{\link{igraph}} undirected graph \emph{g}  or its adjacency matrix (can be extracted from \emph{g}   using \code{\link{get.adjacency}}). Since 3.6.0, if the matrix is provided, such matrix can contain also NAs and the position of such entries will be preserved by the SA  }
  \item{max.iter}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{(2d^3-6d^2+2d+2)} \ln{(e-de)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}      
  \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
          \item{verbose}{ TRUE (default) boolean value. If TRUE print a processing bar during the rewiring algorithm.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}

} 
\details{Performs at most \eqn{max.iter} number of rewiring steps producing a rewired version of an initial undirected graph.
}
\value{
Adjacency matrix of the rewired graph or the relative \emph{igraph} object depending on the input type.
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr


Iorio, F. and and Bernardo-Faura, M. and Gobbi, A. and Cokelaer, T.and Jurman, G.and Saez-Rodriguez, J. (2016) \emph{Efficient randomization of biologicalnetworks while preserving functionalcharacterization of individual nodes}  Bioinformatics 2016 1 (17):542 doi: 10.1186/s12859-016-1402-1.\cr

 Gobbi, A. and Jurman, G. (2013) \emph{Theoretical and algorithmic solutions for null models in network theory} (Doctoral dissertation) \url{http://eprints-phd.biblio.unitn.it/1125/}  \cr
R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
Special thanks to:Davide Albanese
}
\examples{




library(igraph)
library(BiRewire)
g <- erdos.renyi.game(1000,0.1)
##gets the incidence matrix of g
 m<-as.matrix(get.adjacency(graph=g,sparse=FALSE))

## sets parameters
step=1000
max=100*length(E(g))

##rewiring 
m2=birewire.rewire.undirected(m,100*length(E(g))) 
##creates the corresponding bipartite graph 
g2<-graph.adjacency(m2,mode="undirected")


}
\keyword{undirected graph, rewire}
