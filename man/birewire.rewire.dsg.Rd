\encoding{UTF-8}
\name{birewire.rewire.dsg}

\alias{birewire.rewire.dsg}
\title{Efficient rewiring of directed signed graphs}
\description{Optimal implementation of the switching algorithm. It returns the rewired version of the initial directed signed graph.} 
\usage{birewire.rewire.dsg(dsg, max.iter.pos="n", max.iter.neg="n",accuracy=1,verbose=TRUE,
MAXITER_MUL=10,exact=T,save.file=F,path=NA)}
\arguments{
 \item{dsg}{A dsg object: is a list of two incidence matrices (see References), "positive" and "negative", encoding the positive edges and negative edges. This list can be obtained reading a SIF file using \code{getbirewire.load.dsg} function and converting the resulting dataframe using \code{birewire.induced.bipartite}.}
  \item{max.iter.pos}{ "n" (default)  the number of switching steps to be performed on the positive part of  \emph{dsg} (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{(e-de)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{(e-de)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}   
  \item{max.iter.neg}{ "n" (default)  the number of switching steps to be performed on the negative part of  \emph{dsg} (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{(e-de)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{(e-de)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}   
    \item{accuracy}{ 1 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point;}
    \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} \strong{successful} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
  \item{save.file}{ FALSE (default). If TRUE the dsg is converted in a SIF format file and saved in the \emph{path} directory;}
   \item{path}{ NA (default). If \emph{save.file} is true, the dsg is saved in \emph{path} in SIF format. }
} 
\details{ This fuction runs \code{birewire.rewire.bipartite} on the positive and negative part of \emph{dsg}. See references for more details.
}
\value{
Rewired dsg.}
\references{
 Iorio, F. and Gobbi, A.  and  Saez-Rodriguez, J. (2015) \emph{XXX} YYY.\cr

Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi:  <gobbi.andrea@mail.com>\cr
Special thanks to:
Davide Albanese
}
\examples{


library(BiRewire)
data(test_dsg)
dsg=birewire.induced.bipartite(test_dsg)
tmp= birewire.rewire.dsg(dsg,verbose=F)

}
\keyword{directed graph, rewire, pathway, signaling}
