#include <stdlib.h>
#include <R_ext/Print.h>
#include <R.h>
double similarity_undirected(short **m,short **n,size_t ncol,size_t nrow,size_t e);
double similarity(short **m,short **n,size_t ncol,size_t nrow,size_t e );
size_t analysis_ex(short **incidence,size_t ncol, size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose,size_t MAXITER);
size_t rewire_bipartite_ex(short **matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose,size_t MAXITER);
size_t rewire_sparse_bipartite_ex(size_t *from,size_t *to,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose,size_t MAXITER);
size_t analysis_undirected_ex(short **incidence,size_t ncol, size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose,size_t MAXITER);
size_t rewire_ex(short **matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose,size_t MAXITER);
size_t rewire_sparse_ex(size_t *from,size_t *to,size_t *degree,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose,size_t MAXITER);

size_t analysis(short **incidence,size_t ncol, size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose);
size_t rewire_bipartite(short **matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose);
size_t rewire_sparse_bipartite(size_t *from,size_t *to,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose);
size_t analysis_undirected(short **incidence,size_t ncol, size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose);
size_t rewire(short **matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose);
size_t rewire_sparse(size_t *from,size_t *to,size_t *degree,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose);


