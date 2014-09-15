
\encoding{UTF-8}
\name{birewire.analysis}

\alias{birewire.analysis}

\title{Analysis of Jaccard similarity trends across switching steps.}
\description{This function performs a sequence of \emph{max.iter} switching steps on the input bipartite graph \emph{g} and compute the Jaccard similarity between \emph{g} (the initial network) and its rewired version each \emph{step} switching steps. }
\usage{birewire.analysis(incidence, step=10, max.iter="n",accuracy=1,
verbose=TRUE,MAXITER_MUL=10,exact=F)}
\arguments{
  \item{incidence}{Incidence matrix of the initial bipartite graph \emph{g} (can be extracted from an {\code{igraph}} bipartite graph using the {\code{get.incidence}})} function;
  \item{step}{  10 (default): the interval (in terms of switching steps) at which the Jaccard index between \emph{g} and the its current rewired version is computed;}
  \item{max.iter}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of successful switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{(e-de)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{(e-de)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}
 
    \item{accuracy}{ 1 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point;}
    \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation;}
    \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}
 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} successful swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
} 
\details{
 This function performs  \emph{max.iter} switching steps (see references). In particular, at each step two edges are randomly selected from the current version of \emph{g}. Let these two edges be \eqn{(a,b)} and \eqn{(c,d)} (where \eqn{a} and \eqn{c} belong to the first class of nodes whereas \eqn{b} and \eqn{d} belong to the second one), with \eqn{a\not=c} and  \eqn{b\not=d}.\cr 
 If the \eqn{(a,d)} and \eqn{(c,b)} edges are not already present in the current current version of
 \emph{g} then \eqn{(a,d) } and\eqn{ (c,b)} replace \eqn{(a,b)} and \eqn{ (c,d)}. \cr
\cr
At each \emph{step} number of switching steps the function computes the \strong{Jaccard index} between the original graph \emph{g} and its current version.  \cr

}
\value{
A list containing a vector of Jaccard index values computed each (\emph{scores}) switching steps, whose length is equal to \eqn{max.iter/step}, and the analytically derived lower bound (\eqn{N}) of switching steps to be performed by the switching algorithm in order to provide the revired version of \emph{g} with the maximal level of achievable randomness (in terms of dissimilarity from the initial \emph{g}).
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr

  David J. Rogers and Taffee T. Tanimoto, \emph{A Computer Program for Classifying Plants}, Science Vol 132 pp 1115-1118, October 1960 \cr

Hamming, Richard W. (1950), \emph{Error detecting and error correcting codes}, Bell System Technical Journal 29 (2): 147–160, MR 0035935.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
Special thanks to:\cr
Davide Albanese
}
\examples{

library(igraph)
library(BiRewire)
g <- simplify(graph.bipartite( rep(0:1,length=100),
c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
seq(1,100,17),seq(1,100,19),seq(1,100,23),100
)))

##get the incidence matrix of g
 m<-as.matrix(get.incidence(graph=g))

## set parameters
step=1
max=100*length(E(g))

## perform two different analysis using two different maximal number of switching steps
scores<-birewire.analysis(m,step,max)
scores2<-birewire.analysis(m,step,"n")

## plot the Jaccard index scores across intervals of switching steps
plot(x=step*seq(1:length(scores$similarity_scores)),y= scores$similarity_scores,
type='l',xlab="Number of rewiring",ylab="Index value",ylim=c(0,1))

lines(x=step*seq(1:length(scores2$similarity_scores)),y= scores2$similarity_scores,
col="red")
legend(max*0.5,1, c("Jaccard index","Jaccard index with lower-bound N"), cex=0.9, 
col=c("black","red"), lty=1:1,lwd=3)


}
\keyword{rewire, bipartite graph}
