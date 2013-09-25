\name{BRCA_binary_matrix}
\docType{data}
\alias{BRCA_binary_matrix}
\title{TCGA Brest Cancer data}
\description{
Breast cancer samples and their respective mutations downloaded from the Cancer Cancer Genome Atlas (TCGA), used in \cite{Gobbi er al.}.
Germline mutations were filtered out of the list of reported mutations; synonymous mutations and mutations identified as benign and tolerated were also removed from the dataset. The bipartite graph resulting when considering this matrix as an incidence matrix has \eqn{n_r=757, n_c=9757, e=19758}  for an edge density equal to \eqn{0.27\%}.

}
\references{
 Gobbi, A. and Iorio, and Wedge, D. and Dawson, K. and Ludmil, A. F. and Jurman, G. and Saez-Rodriguez, J.  (submitted), \emph{Preserving genomic-event distributions in
randomized next generation sequencing data}.\cr


  }
\usage{data(BRCA_binary_matrix)}
\source{http://tcga.cancer.gov/dataportal/}
\keyword{datasets}
