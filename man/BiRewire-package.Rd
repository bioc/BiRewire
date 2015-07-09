\encoding{UTF-8}
\name{BiRewire-package}
\alias{BiRewire-package}
\alias{BiRewire}
\docType{package}
\title{The BiRewire package}
\description{
   R package for computationally-efficient rewiring of bipartite graphs (or randomisation of 0-1 tables with prescribed marginal totals) and undirected. The package provides useful functions for the analysis and the randomisation of large biological datasets that can be encoded as 0-1 tables, hence modeled as bipartite graphs by considering a 0-1 table as an incidence matrix. Large collections of such randomised tables can be used to approximate null models, preserving event-rates both across rows and columns, for statistical significance tests of combinatorial properties of the original dataset. The package provides an interface to a sampler routine useful for generating correctly such collections. Moreover a visual monitoring for the Markov Chain underlying the swithicng algorithm has been implemented.

   }
\details{
Summary:
  \tabular{ll}{
  Package: \tab BiRewire\cr
  Version: \tab 2.2.0\cr
  Date: 2015-07-07\cr
  Require: slam, igraph,tsne,  R>=2.10\cr
  URL:  http://www.ebi.ac.uk/~iorio/BiRewire  \cr
  License: \tab GPL-3 \cr
  }
  
  
 
}
\author{
Andrea Gobbi [aut], Davide Albanese [cbt], Francesco Iorio [cbt], Giuseppe Jurman [cbt].

  
  Maintainer: Andrea Gobbi <gobbi.andrea@mail.com>

}
\references{

Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr
Jaccard, P. (1901), \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura},
 Bulletin de la Société Vaudoise des Sciences Naturelles 37: 547–579.\cr
\cr
R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon (2003), \emph{On the uniform generation of random graphs with prescribed degree sequences}, eprint arXiv:cond-mat/0312028 
 Csardi, G. and  Nepusz, T (2006)\cr
  Van der Maaten, L.J.P. and Hinton, G.E., \emph{ Visualizing High-Dimensional Data Using t-SNE}. Journal of Machine Learning Research 9(Nov):2579-2605, 2008
\emph{The igraph software package for complex network research},
 InterJournal, Complex Systems 
    \url{http://igraph.sf.net}


}
\keyword{ package }
