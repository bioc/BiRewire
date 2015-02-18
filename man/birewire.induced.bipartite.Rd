\encoding{UTF-8}
\name{birewire.induced.bipartite}

\alias{birewire.induced.bipartite}
\title{Transform sif data frame into positive and negarive incidnece matrix.}
\description{The routine transforms the initial dsg graph (in SIF) format into two bipartite graphs: one for positive edges and the other for negative edges.} 
\usage{birewire.induced.bipartite(g,delimitators=list(negative='-',positive='+'))}
\arguments{
 \item{g}{A dataframe in SIF format describing a dsg (for example the output of \code{birewire.load.dsg};}
  \item{delimitators}{list(negative='-',positive='+') (default):a list with 'positive' and 'negative' names identifying the character encoding the relation;}


 } 
\details{ This fuction extract the positive and negative part of \emph{g} and create a dsg object that can be used for example in the rewiring algorithm. See references for more details. Is is the inverse function of \code{birewire.build.dsg}.
}
\value{
A list of two incidence matrix.}
\references{
 Iorio, F. and Gobbi, A.  and  Saez-Rodriguez, J. (2015) \emph{XXX} YYY.\cr
  }

\examples{

data(test_dsg)
dsg=birewire.induced.bipartite(test_dsg)

}
\keyword{directed graph, rewire, pathway, signaling}
