\encoding{UTF-8}
\name{birewire.rewire}

\alias{birewire.rewire}
\title{Efficient rewiring of undirected graphs}
\description{Optimal implementation of the switching algorithm. It returns the rewired version of the initial undirected graph or its adjacency matrix.} 
\usage{birewire.rewire(adjacency, max.iter="n",accuracy=1,
verbose=TRUE,MAXITER_MUL=10,exact=F)}
\arguments{
  \item{adjacency}{An \code{igraph} undirected graph \emph{g}  or its adjacency matrix (can be extracted from \emph{g}   using {\code{get.adjacency(g)}});  }
  \item{max.iter}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{(2d^3-6d^2+2d+2)} \ln{(e-de)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{(e-de)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}      
   \item{accuracy}{ 1 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point;}
        \item{verbose}{ TRUE (default) boolean value. If TRUE print a processing bar during the rewiring algorithm.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \code{max.iter} \strong{successful} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
} 
\details{Performs at most \eqn{max.iter} number of rewiring steps producing a rewired version of an initial undirected graph.
}
\value{
Adjacency matrix of the rewired graph. 
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr
 Gobbi, A. and Jurman, G.  (in preparation),  \emph{Number of required Switching Steps in the Switching Algorithm for undirected graphs}.\cr
R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi
\cr
Special thanks to:\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
Davide Albanese
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
m2=birewire.rewire(m,100*length(E(g))) 
##creates the corresponding bipartite graph 
g2<-graph.adjacency(m2,mode="undirected")


}
\keyword{undirected graph, rewire}
