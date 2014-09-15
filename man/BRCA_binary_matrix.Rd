\name{BRCA_binary_matrix}
\docType{data}
\alias{BRCA_binary_matrix}
\title{TCGA Brest Cancer data}
\description{
Breast cancer samples and their respective mutations downloaded from the Cancer Cancer Genome Atlas (TCGA), used in \cite{Gobbi et al.}.
Germline mutations were filtered out of the list of reported mutations; synonymous mutations and mutations identified as benign and tolerated were also removed from the dataset. The bipartite graph resulting when considering this matrix as an incidence matrix has \eqn{n_r=757, n_c=9757, e=19758}  for an edge density equal to \eqn{0.27\%}.

}
\references{
Gobbi, A. and  Iorio, F. and  Dawson, K. J. and Wedge, D. C. and Tamborero, D. and Alexandrov, L. B. and Lopez-Bigas, N. and  Garnett, M. J. and  Jurman, G. and  Saez-Rodriguez, J. (2014) \emph{Fast randomization of large genomic datasets while preserving alteration counts} Bioinformatics 2014 30 (17): i617-i623 doi: 10.1093/bioinformatics/btu474.\cr


  }
\usage{data(BRCA_binary_matrix)}
\source{http://tcga.cancer.gov/dataportal/}
\keyword{datasets}
