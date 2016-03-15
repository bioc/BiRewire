\encoding{UTF-8}
\name{birewire.rewire.bipartite}

\alias{birewire.rewire.bipartite}
\title{Efficient rewiring of bipartite graphs}
\description{Optimal implementation of the switching algorithm. It returns the rewired version of the initial bipartite graph or its incidence matrix.} 
\usage{birewire.rewire.bipartite(incidence, max.iter="n",accuracy=0.00005,verbose=TRUE,
MAXITER_MUL=10,exact=FALSE)}
\arguments{
 \item{incidence}{Incidence matrix of the initial bipartite graph \emph{g} (can be extracted from an \code{\link{igraph}} bipartite graph using the \code{\link{get.incidence}}) function; or the entire bipartite \code{\link{igraph}} graph}
  \item{max.iter}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{((e-de)/\delta)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}   

 \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
         \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}

} 
\details{Main function of the package. It performs at most \eqn{max.iter} switching steps producing a rewired version of an initial bipartite graph.
}
\value{
Incidence matrix  of the rewired graphn or the \emph{igraph} corresponding object depending on the input type.
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr

}
\examples{



library(igraph)
library(BiRewire)
g <-graph.bipartite( rep(0:1,length=10), c(1:10))

##gets the incidence matrix of g
 m<-as.matrix(get.incidence(graph=g))

##rewiring
m2=birewire.rewire.bipartite(m,100*length(E(g))) 
##creates the corresponding bipartite graph 
g2<-birewire.bipartite.from.incidence(m2,directed=TRUE)


}
\keyword{bipartite graph, rewire}
