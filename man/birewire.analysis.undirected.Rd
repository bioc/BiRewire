\encoding{UTF-8}
\name{birewire.analysis.undirected}

\alias{birewire.analysis.undirected}

\title{Analysis of Jaccard similarity trends across switching steps.}
\description{This function performs a sequence of \emph{max.iter} switching steps on the input undirected graph \emph{g} and compute the Jaccard similarity between \emph{g} (the initial network) and its rewired version each \emph{step} switching steps. This procedure is pefromed \emph{n.networks} times and a simple explorative plot, with mean and CI, is visualized if \emph{display} is set to true. }
\usage{birewire.analysis.undirected(adjacency, step=10, max.iter="n",accuracy=0.00005,
verbose=TRUE,MAXITER_MUL=10,exact=FALSE,n.networks=50,display=TRUE)}
\arguments{
    \item{adjacency}{Incidence matrix of the initial bipartite graph \emph{g} (can be extracted from an {\code{igraph}} undirected graph using the \code{get.adjacency} function)};
  \item{step}{  10 (default): the interval (in terms of switching steps) at which the Jaccard index between \emph{g} and the its current rewired version is computed;}
  \item{max.iter}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of successful switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{(2d^3-dd^2+2d+2)} \ln{((e-de)/\delta)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}
 
 \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
         \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation;}
    \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}
  \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
 \item{n.networks}{ 50 (default), the number of independent rewiring process starting from the same inital graph from which the mean value and the CI is computed. }
  \item{display}{ TRUE (default). If TRUE two explorative plots are displayed summarizing the trend of the Jaccard index in terms of mean and confidence interval.}
} 
\details{
 This function performs  \emph{max.iter} switching steps (see references). In particular, at each step two edges are randomly selected from the current version of
 \emph{g}. Let these two edges be \eqn{(a,b)} and \eqn{(c,d)}, with \eqn{a\not=c}, \eqn{b\not=d},  \eqn{a\not=d}, \eqn{b\not=c} .\cr 
 If the \eqn{(a,d)} and \eqn{(c,b)} (or  \eqn{(a,d)} and \eqn{(b,d)}) edges are not already present in the current version of
 \emph{g} then \eqn{(a,d) } and \eqn{ (c,b)} replace \eqn{(a,b)} and \eqn{ (c,d)} (or \eqn{(a,b)} and \eqn{ (c,d)} replace \eqn{(a,c) } and \eqn{ (b,d)}). If both of the configuarations are allowed, then one of them is randomly selected. \cr
\cr
At each \emph{step} number of switching steps the function computes the \strong{Jaccard index} between the original graph \emph{g} and its current version.  \cr
This procedure is perfomed \emph{n.networks} times and if \emph{display} is set to TRUE, two explorative plots showing the mean value of the Jaccad Index over the SS and its CI are displayed.

}
\value{
A list containing a data.frame \emph{data} collecting all the Jacard index computed (each row is a run of the SA), and the analytically derived lower bound \emph{N} of switching steps to be performed by the switching algorithm in order to provide the revired version of \emph{g} with the maximal level of achievable randomness (in terms of dissimilarity from the initial \emph{g}).

}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

 Gobbi, A. and Jurman, G.  (in preparation),  \emph{Number of required Switching Steps in the Switching Algorithm for undirected graphs}.\cr

Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
Special thanks to:\cr
Davide Albanese
}
\examples{

library(BiRewire)
g <- erdos.renyi.game(1000,0.1)
##get the incidence matrix of g
 m<-as.matrix(get.adjacency(graph=g,sparse=FALSE))

## set parameters
step=1000
max=100*length(E(g))

## perform two different analysis using two different numbers of switching steps
scores<-birewire.analysis.undirected(m,step,max,n.networks=10)
scores2<-birewire.analysis.undirected(m,step,"n",n.networks=10)


}
\keyword{rewire, undirected graph}
