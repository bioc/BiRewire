\encoding{UTF-8}
\name{birewire.visual.monitoring.undirected}

\alias{birewire.visual.monitoring.undirected}

\title{Visual monitoring of the Markov chain underlying the SA for undirected graphs.}
\description{This function generates a cascade-sampling from the model at different switching steps given in \emph{sequence}. For each step the routine computes the pairwise Jaccard distance (1-JI) among the samples and perfroms, on the resulting matix, a dimentional scaling reduction (using \code{\link{Rtsne}}).  If \emph{display} is set to \emph{TRUE} the relative plot is displayed.  }
\usage{birewire.visual.monitoring.undirected(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,
	exact=FALSE,n.networks=100,perplexity=15,sequence=c(1,5,100,"n"),ncol=2,
	nrow=length(sequence)/ncol,display=TRUE)}
\arguments{
  \item{data}{The initial undirected graph, either an adjacency matrix or an \emph{igraph} undirected graph object.  Since 3.6.0, if the matrix is provided, such matrix can contain also NAs and the position of such entries will be preserved by the SA;}
  \item{accuracy}{ 0.00005 (default) is the desired level of accuracy reflecting the average distance between the Jaccard index at the N-th step and its analytically derived fixed point in terms of fracion of common edges;}
         \item{verbose}{ TRUE (default). When TRUE a progression bar is printed during computation.}
 \item{MAXITER_MUL}{ 10 (default). If \emph{exact==TRUE}  in order to prevent a possible infinite loop the program stops anyway after MAXITER_MUL*max.iter iterations;}

 \item{exact}{ FALSE (default). If TRUE the program performs \emph{max.iter} swithcing steps, otherwise the program will count also the not-performed swithcing steps;}

  \item{n.networks}{  100 (default): the number of network generated for each step defined in \emph{sequence} ;}
  \item{perplexity}{ 15 (default): the value of perplexity passed to the function \code{\link{Rtsne}};}
  \item{sequence}{ c(1,5,100,"n")(default) the sequence of step for wich generating a sampler 
(see \code{\link{birewire.sampler.undirected}})}
  \item{ncol}{ 2 (default). The number of column in the plot;}
  \item{nrow}{ length(sequence)/ncol (default). The number of row in the plot;}
  \item{display}{ TRUE (default). If TRUE the result of tsne is displayed.}
} 
\details{
 For each value \emph{p} in \emph{sequence} (it that can also contain the special character "n", see \code{\link{birewire.rewire.bipartite}}), the routine generates \emph{n.networks} sampled each \emph{p} SS from the SA initialized with the given \emph{data}. Pariwise distance are computed using the Jaccard distance and the resulting matrix is the input for the dimensional scaling performed by the function \code{\link{Rtsne}}. An explorative plot is displayed if \emph{display} is set to TRUE. 
}
\value{
A list containing the list containing the distance matrices  \emph{dist} and the list containing the tsne results \emph{Rtsne}.
}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr

Iorio, F. and and Bernardo-Faura, M. and Gobbi, A. and Cokelaer, T.and Jurman, G.and Saez-Rodriguez, J. (2016) \emph{Efficient randomization of biologicalnetworks while preserving functionalcharacterization of individual nodes}  Bioinformatics 2016 1 (17):542 doi: 10.1186/s12859-016-1402-1.\cr

Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr

R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028

Van der Maaten, L.J.P. and Hinton, G.E. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008
  }
\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>\cr
}
\examples{
library(BiRewire)
g <- erdos.renyi.game(1000,0.1)
tsne = birewire.visual.monitoring.undirected(g,display=FALSE,n.networks=10,perplexity=1)


}
\keyword{rewire, bipartite graph}
