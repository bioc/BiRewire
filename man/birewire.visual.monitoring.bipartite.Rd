\encoding{UTF-8}
\name{birewire.visual.monitoring.bipartite}

\alias{birewire.visual.monitoring.bipartite}

\title{Visual monitoring of the Markov chain underlying the SA for directed graphs.}
\description{This function generates a cascade-sampling from the model at different switching steps given in \emph{sequence}. For each step the routine computes the pairwise Jaccard distance (1-JI) among the samples and perfroms, on the resulting matix, a dimentional scaling reduction (using \code{tsne}). If \emph{display} is true, the relative plot is displayed.  }
\usage{birewire.visual.monitoring.bipartite(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,
  exact=FALSE,n.networks=100,perplexity=15,sequence=c(1,5,100,"n"),ncol=2,
  nrow=length(sequence)/ncol,display=TRUE)}
\arguments{
  \item{data}{The initial bipartite graph, either an incidence matrix or an \emph{igraph} bipartite graph object;}
\item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
         \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}
  \item{n.networks}{  100 (default): the number of network generated for each step defined in \emph{sequence} ;}
  \item{perplexity}{ 15 (default): the value of perplexity passed to the function \code{tsne};}
  \item{sequence}{ c(1,5,100,"n")(default) the sequence of step for wich generating a sampler (see\code{birewire.rewire.sampler.undirected})}
  \item{ncol}{ 2 (default). The number of column in the plot;}
  \item{nrow}{ length(sequence)/ncol (default). The number of row in the plot;}
  \item{display}{ TRUE (default). If TRUE the result is displayed.} 
  }
\details{
 For each value \emph{p} in \emph{sequence} (it that can also contain the special character "n", see \code{birewire.rewire.bipartite}), the routine generates \emph{n.networks} sampled each \emph{p} SS from the SA initialized with the given \emph{data}. Pariwise distance are computed using the Jaccard distance and the resulting matrix is the input for the dimensional scaling performed by the function \code{tsne}. An explorative plot is displayed if \emph{display} is set to TRUE. 
}
\value{
A list containing the list containing the distance matrices  \emph{dist} and the list containing the tsne results \emph{tsne}.
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028

Van der Maaten, L.J.P. and Hinton, G.E. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
Special thanks to:\cr
Davide Albanese
}
\examples{

library(BiRewire)
g <- graph.bipartite( rep(0:1,length=100), c(1:100))
birewire.visual.monitoring.bipartite(g,display=FALSE,n.networks=10)


}
\keyword{rewire, bipartite graph}
