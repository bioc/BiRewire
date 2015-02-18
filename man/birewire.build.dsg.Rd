\encoding{UTF-8}
\name{birewire.build.dsg}

\alias{birewire.build.dsg}
\title{Transform a dsg  into SIF format.}
\description{The routine transforms the initial dsg obgect (two bipartite graphs) into SIF dsg format.} 
\usage{birewire.build.dsg(dsg,delimitators=list(negative='-',positive='+'))}
\arguments{
 \item{dsg}{The dsg to be converted;}
  \item{delimitators}{list(negative='-',positive='+') (default):a list with 'positive' and 'negative' names identifying the character encoding the relation;}

 } 
\details{ This fuction converts the dsg object into a SIF format that can be saved using \code{bipartite.write.dsg} using the given delimitators for encoding the relations. See references for more details. Is is the inverse function of \code{birewire.induced.bipartite}.
}
\value{
A dsg in SIF format.}
\references{
 Iorio, F. and Gobbi, A.  and  Saez-Rodriguez, J. (2015) \emph{XXX} YYY.\cr
  }

\examples{

data(test_dsg)
dsg=birewire.induced.bipartite(test_dsg)
tmp= birewire.rewire.dsg(dsg,verbose=FALSE)
dsg2=birewire.build.dsg(tmp)

}
\keyword{directed graph, rewire, pathway, signaling}
