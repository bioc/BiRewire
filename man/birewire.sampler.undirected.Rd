\encoding{UTF-8}
\name{birewire.sampler.undirected}

\alias{birewire.sampler.undirected}
\title{Efficient generation of a null model for a given undirected graph}
\description{ The routine samples correctly from the null model of a given undirected graph creating a set of randomized version of the initial undirected graph.} 
\usage{birewire.sampler.undirected(adjacency,K,path,max.iter="n", accuracy=0.00005,
	verbose=TRUE,MAXITER_MUL=10,exact=FALSE,write.sparse=TRUE)

}
\arguments{
 \item{adjacency}{Adjacency matrix of the initial undirected graph;}
 \item{K}{  The number of networks that has to be generated;}
   \item{path}{ The directory in which the routine stores the outputs; }

  \item{max.iter}{ "n" (default)  see \code{\link{birewire.rewire.undirected}} for references}
  \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}  
    \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
 \item{write.sparse}{ TRUE (default). If FALSE the table is written as an R data.frame (long time and more space needed)}

} 
\details{ The routine creates, starting from the given path, different subfolders in order to have maximum 1000 files for folder . Moreover the incidence matrices are saved using \code{\link{write_stm_CLUTO}} (sparse matrices) that can be loaded using \code{\link{read_stm_CLUTO}}. The set is generated calling birewire.rewire.undirected on the last generated graph starting from the input graph.
}	
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

Gobbi, A. and Jurman, G. (2013) \emph{Theoretical and algorithmic solutions for null models in network theory} (Doctoral dissertation) \url{http://eprints-phd.biblio.unitn.it/1125/}  \cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028
  }
\author{Andrea Gobbi:  <gobbi.andrea@mail.com>\cr

}

\keyword{undirected graph, rewire	}
