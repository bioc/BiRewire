\encoding{UTF-8}
\name{birewire.analysis.dsg}

\alias{birewire.analysis.dsg}

\title{Analysis of Jaccard similarity trends across switching steps.}
\description{This function performs a sequence of \emph{max.iter.pos} (and \emph{max.iter.pos}) switching steps on the positive (and negative) part of the input dsg \emph{g} and computes the Jaccard similarity between \emph{g} (the initial network) and its rewired version each \emph{step} switching steps. This procedure is pefromed \emph{n.networks} times and a simple explorative plot, with mean and CI, is visualized if \emph{display} is set to \emph{TRUE}. The plot shows the trend of the Jaccad Index relative to the positve (and negative) part of \emph{g}. }
\usage{birewire.analysis.dsg(dsg, step=10, max.iter.pos='n',max.iter.neg='n',accuracy=0.00005,
			     verbose=TRUE,MAXITER_MUL=10,exact=FALSE,n.networks=50,display=TRUE)}
\arguments{
    \item{dsg}{The initial dsg object (see \code{\link{birewire.induced.bipartite}}). Note that the dsg must contain a list of two incidence matrices and not \code{\link{igraph}} bipartite graphs. }
  \item{step}{  10 (default): the interval (in terms of switching steps) at which the Jaccard index between \emph{g} and the its current rewired version is computed;}
  \item{max.iter.pos}{ "n" (default)  the number of switching steps to be performed (or if \emph{exact==TRUE} the number of successful switching steps) for the positive part of \emph{g}.
 See \code{\link{birewire.rewire.bipartite}} for more details;}
   \item{max.iter.neg}{ "n" (default)  the same of \emph{max.iter.p} but relative to the negative part;}
 \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
         \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation;}
    \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}
  \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
 \item{n.networks}{ 50 (default), the number of independent rewiring process starting from the same inital graph from which the mean value and the CI is computed. }
  \item{display}{ TRUE (default). If TRUE two explorative plots are displayed summarizing the trend of the Jaccard index in terms of mean and confidence interval.}
} 
\details{
 This procedure acts in the same way of \code{\link{birewire.analysis.bipartite}} but in the case of dsg. The similarity is measurwe using \code{\link{birewire.similarity.dsg}}.

 }
\value{
A list containing two lists: \emph{data} that is a list collecting all the Jacard index computed (each row is a run of the SA) for the positive and negative part, and a list with the analytically derived lower bounds \emph{N} for the positive and negative part of \emph{g}.
}
\references{

 Iorio, F. and Gobbi, A.  (2016) \emph{In Preparation} In Preparation.\cr

Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
}
\examples{

library(BiRewire)
data(test_dsg)
dsg <-  birewire.induced.bipartite(test_dsg,sparse=FALSE)

a=birewire.analysis.dsg(dsg,verbose=FALSE,step=1,exact=TRUE,max.iter.pos=200,max.iter.neg=50)

}
\keyword{rewire,dsg}
