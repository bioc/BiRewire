\encoding{UTF-8}
\name{birewire.rewire.bipartite.and.projections}

\alias{birewire.rewire.bipartite.and.projections}

\title{Analysis and rewiring function processing a bipartite graphs and its two projections}
\description{This function performs the same analysis of \code{\link{birewire.analysis.bipartite}} but additionally it provides in output a rewired version of the two networks resulting from the natural projections of the initial graph, together with the corresponding Jaccard index trends.} 
\usage{birewire.rewire.bipartite.and.projections(graph,step=10,max.iter="n",
accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10)}
\arguments{
  \item{graph}{ A bipartite graph \emph{g};}
    \item{max.iter}{ "n" (default)  the  number of successful switching steps to be performed. If equal to "n" then this number is considered equal to the analytically derived lower bound \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } presented in \emph{Gobbi et al.} (see References); }
  \item{step}{  10 (default): the interval (in terms of switching steps) at which the Jaccard index between \emph{g} and the its current rewired version is computed;}
 \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
             \item{verbose}{ TRUE (default) boolean value. If TRUE print a processing bar during the rewiring algorithm.}
   \item{MAXITER_MUL}{ 10 (default).Since \eqn{N} indicates the number of successful switching steps, in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations  ;}
} 
\details{
See \code{\link{birewire.analysis.bipartite}} for details.
}
\value{
A list containing the three sequences of Jaccard index values (similarity_scores, similarity_scores.proj1, similarity_scores.proj2) for the three resulting graphs respectively (rewired, rewired.proj1, rewired.proj2). The first one is the rewired
version of the initial graph \emph{g}, while the second and the third one are rewired versions of its natural projections.
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

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
seq(1,100,17),seq(1,100,19),seq(1,100,23),100)))
##gets the incidence matrix of g
 m<-as.matrix(get.incidence(graph=g))
## rewires g and its projections
result=birewire.rewire.bipartite.and.projections(g,step=10,max.iter="n",accuracy=0.00005) 



}

\keyword{bipartite graph, projection, rewire}
