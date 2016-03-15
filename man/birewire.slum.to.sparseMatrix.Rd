\encoding{UTF-8}
\name{birewire.slum.to.sparseMatrix}

\alias{birewire.slum.to.sparseMatrix}

\title{The function transforms a triplet sparse matrix from slum package to a Matrix sparse matrix.}
\description{Transform a triplet sparse matrix from \emph{slum} package to a \emph{Matrix} sparse matrix that can be used by \emph{igraph} for creating a network. This function could be used in order to analyze graphs obtained from samplers routines (\code{\link{birewire.sampler.undirected}},\code{\link{birewire.sampler.dsg}} and \code{\link{birewire.sampler.bipartite}}.)}  
\usage{birewire.slum.to.sparseMatrix( simple_triplet_matrix_sparse)}
\arguments{
  \item{simple_triplet_matrix_sparse}{A triplet sparse matrix, usually the object coming from \code{\link{read_stm_CLUTO}}.}
    
  }
\value{
Returns an Matrix sparse matrix that could be used for building an \code{\link{igraph}} graph using \code{\link{graph.adjacency}}.
}

\author{Andrea Gobbi
\cr
Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>
}


\keyword{slum, Matrix,sparse matrix}
