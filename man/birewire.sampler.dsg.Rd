\encoding{UTF-8}
\name{birewire.sampler.dsg}

\alias{birewire.sampler.dsg}
\title{Efficient generation of a null model for a given dsg.}
\description{Efficient generation of a null model for a given dsg. The routine samples correctly from the null model of a given dsg creating a set of randomized dsgs.} 
\usage{birewire.sampler.dsg(dsg,K,path,delimitators=list(negative='-',positive='+'),exact=FALSE,
  verbose=TRUE, max.iter.pos='n',max.iter.neg='n', accuracy=0.00005,MAXITER_MUL=10,check_pos_neg=FALSE)
}
\arguments{
 \item{dsg}{A dsg object: is a list of two incidence matrices (see References), "positive" and "negative", encoding the positive edges and negative edges. This list can be obtained reading a SIF file using \code{\link{birewire.load.dsg}} function and converting the resulting dataframe using \code{\link{birewire.induced.bipartite}}.}
  \item{max.iter.pos}{ "n" (default)  the number of switching steps to be performed on the positive part of  \emph{dsg} (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{((e-de)/\delta)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}   
  \item{max.iter.neg}{ "n" (default)  the number of switching steps to be performed on the negative part of  \emph{dsg} (or if \emph{exact==TRUE} the number of \strong{successful} switching steps). If equal to "n" then this number is considered equal to the analytically derived lower bound presented in \emph{Gobbi et al.} (see References): \eqn{N={e}/{2(1-d)} \ln{((e-de)/\delta)} } if exact is FALSE, \eqn{N={e(1-d)}/{2} \ln{((e-de)/\delta)} } otherwise ,  
  where \eqn{e}  is the number of edges of \emph{g} and \eqn{d} its edge density . This bound is much lower than the empirical one proposed  in  \emph{Milo et al. 2003} (see References);}   
    \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
    \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}

    \item{path}{ The directory in which the routine stores the outputs; }
      \item{K}{  The number of network that has to be generated;}
  \item{delimitators}{list(negative='-',positive='+') (default):a list with 'positive' and 'negative' names identifying the character encoding the relation used for writing the ouput with \code{\link{birewire.build.dsg}};}
    \item{check_pos_neg}{FALSE (default) : if true check if there are positive and negative loops in the same arc, if so the routine will not save the dsg and goes on with the marcov chain (it can be slower due to the chekcs);}

} 
\details{ The routine creates, starting from a given path, different subfolders in order to have maximum 1000 files for folder; the SIF files are saved using \code{birewire.write.dsg}, an internal routine. The set is generated calling birewire.rewire.dsg on the last generated dsg starting from the input dsg.
}
\references{

Iorio, F. and and Bernardo-Faura, M. and Gobbi, A. and Cokelaer, T.and Jurman, G.and Saez-Rodriguez, J. (2016) \emph{Efficient randomization of biologicalnetworks while preserving functionalcharacterization of individual nodes}  Bioinformatics 2016 1 (17):542 doi: 10.1186/s12859-016-1402-1.\cr

Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi:  <gobbi.andrea@mail.com>\cr

}

\keyword{directed graph, rewire, pathway, signaling}


