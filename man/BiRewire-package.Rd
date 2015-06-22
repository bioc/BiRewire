\encoding{UTF-8}
\name{BiRewire-package}
\alias{BiRewire-package}
\alias{BiRewire}
\docType{package}
\title{The BiRewire package}
\description{
   R package for computationally-efficient rewiring of bipartite graphs (or randomisation of 0-1 tables with prescribed marginal totals) and directed signed graphs (dsg). The package provides useful functions for the analysis and the randomisation of large biological datasets that can be encoded as 0-1 tables, hence modeled as bipartite graphs by considering a 0-1 table as an incidence matrix, and for data that can be encoded as directed signed graphs such as patways and signaling networks. Large collections of such randomised tables can be used to approximate null models, preserving event-rates both across rows and columns, for statistical significance tests of combinatorial properties of the original dataset. The package provides also some sampler routines useful for generating correctly such collections.
Routines for undirected graphs are also provided.

   }
\details{
Summary:
  \tabular{ll}{
  Package: \tab BiRewire\cr
  Version: \tab 1.12.1\cr
  Date: 2015-03-01\cr
  Require: slam, igraph, R>=2.10\cr
  URL:  http://www.ebi.ac.uk/~iorio/BiRewire  \cr
  License: \tab GPL-3 \cr
  }
  
  
 
}
\author{
Andrea Gobbi [aut], Davide Albanese [cbt], Francesco Iorio [cbt], Giuseppe Jurman [cbt].

  
  Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>

}
\references{

 Iorio, F. and Gobbi, A. (2015) \emph{In Preparation} In Preparation.\cr

Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr
Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr
\cr
  David J. Rogers and Taffee T. Tanimoto, \emph{A Computer Program for Classifying Plants}, Science Vol 132 pp 1115-1118, October 1960 \cr
\cr
Hamming, Richard W. (1950), \emph{Error detecting and error correcting codes}, Bell System Technical Journal 29 (2): 147–160, MR 0035935.\cr
\cr
R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028 \cr
\cr
 Csardi, G. and  Nepusz, T (2006)
\emph{The igraph software package for complex network research},
 InterJournal, Complex Systems 
    \url{http://igraph.sf.net}


}
\keyword{ package }
