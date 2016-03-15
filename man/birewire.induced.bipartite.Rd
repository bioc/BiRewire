\encoding{UTF-8}
\name{birewire.induced.bipartite}

\alias{birewire.induced.bipartite}
\title{Transform a SIF data frame into a dsg object  (a list of positive and negative incidence matrix).}
\description{The routine transforms the initial dsg graph in SIF format into a dsg object made of two bipartite graphs: one for positive edges and the other for negative edges.} 
\usage{birewire.induced.bipartite(g,delimitators=list(negative='-',positive='+'),sparse=FALSE)}
\arguments{
 \item{g}{A dataframe in SIF format describing a dsg (for example the output of \code{\link{birewire.load.dsg}};}
\item{delimitators}{list(negative='-',positive='+') (default):a list with 'positive' and 'negative' names identifying the character encoding the relation;}
\item{sparse}{ FALSE (default): if TRUE the two bipartite graphs are saved as \code{\link{igraph}} bipartite graphs ;}

 } 
\details{ This fuction extract the positive and negative part of \emph{g} and create a dsg object that can be used for example in the rewiring algorithm. Is is the inverse function of \code{\link{birewire.build.dsg}}.
}
\value{
A list of two incidence matrix or bipartite \code{\link{igraph}} objects.}

\references{
	
	 Iorio, F. and Gobbi, A.  (2016) \emph{In Preparation} In Preparation.\cr

}
\examples{

data(test_dsg)
dsg=birewire.induced.bipartite(test_dsg)

}
\keyword{directed graph, rewire, pathway, signaling}
