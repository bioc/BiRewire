
\encoding{UTF-8}
\name{birewire.analysis.bipartite}

\alias{birewire.analysis.bipartite}

\title{Analysis of Jaccard similarity trends across switching steps.}
\description{This function performs a sequence of \emph{max.iter} switching steps on the input bipartite graph \emph{g} and compute the Jaccard similarity between \emph{g} (the initial network) and its rewired version each \emph{step} switching steps. This procedure is pefromed \emph{n.networks} times and a simple explorative plot, with mean and CI, is visualized if \emph{display} is set to true. }
\usage{birewire.analysis.bipartite(incidence, step=10, max.iter="n",accuracy=0.00005,
verbose=TRUE,MAXITER_MUL=10,exact=FALSE,n.networks=50,display=TRUE)}
\arguments{
  \item{incidence}{Incidence matrix of the initial bipartite graph \emph{g} (can be extracted from an \code{\link{igraph}} bipartite graph using the \code{\link{get.incidence}} function)};
  \item{step}{  10 (default): the interval (in terms of switching steps) at which the Jaccard index between \emph{g} and the its current rewired version is computed;}
  \item{max.iter}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of successful switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{((e-de)/\delta)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}
 
 \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
         \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation;}
    \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}
  \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
 \item{n.networks}{ 50 (default), the number of independent rewiring process starting from the same inital graph from which the mean value and the CI is computed. }
  \item{display}{ TRUE (default). If TRUE two explorative plots are displayed summarizing the trend of the Jaccard index in terms of mean and confidence interval.}
} 
\details{
 This function performs  \emph{max.iter} switching steps (see references). In particular, at each step two edges are randomly selected from the current version of \emph{g}. Let these two edges be \eqn{(a,b)} and \eqn{(c,d)} (where \eqn{a} and \eqn{c} belong to the first class of nodes whereas \eqn{b} and \eqn{d} belong to the second one), with \eqn{a\not=c} and  \eqn{b\not=d}.\cr 
 If the \eqn{(a,d)} and \eqn{(c,b)} edges are not already present in the current current version of
 \emph{g} then \eqn{(a,d) } and\eqn{ (c,b)} replace \eqn{(a,b)} and \eqn{ (c,d)}. \cr
\cr
At each \emph{step} number of switching steps the function computes the \strong{Jaccard index} between the original graph \emph{g} and its current version.  \cr
This procedure is perfomed \emph{n.networks} times and if \emph{display} is set to TRUE, two explorative plots showing the mean value of the Jaccad Index over the SS and its CI are displayed.
}
\value{
A list containing a data.frame \emph{data} collecting all the Jacard index computed (each row is a run of the SA), and the analytically derived lower bound \emph{N} of switching steps to be performed by the switching algorithm in order to provide the revired version of \emph{g} with the maximal level of achievable randomness (in terms of dissimilarity from the initial \emph{g}).
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
Special thanks to:\cr Davide Albanese
}
\examples{
library(BiRewire)
g <-graph.bipartite(rep(0:1,length=10), c(1:10))

##get the incidence matrix of g
 m<-as.matrix(get.incidence(graph=g))

## set parameters
step=1
max=100*length(E(g))

## perform two different analysis using two different maximal number of switching steps
scores<-birewire.analysis.bipartite(m,step,max,n.networks=10)
scores2<-birewire.analysis.bipartite(m,step,"n",n.networks=10)



}
\keyword{rewire, bipartite graph}
