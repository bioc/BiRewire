/*
 This code is written by Andrea Gobbi  <gobbi.andrea@mail.com> 2013
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "BiRewire.h"

SEXP R_analysis (SEXP incidence,SEXP ncol, SEXP nrow, SEXP step, SEXP max_iter, SEXP verbose, SEXP MAXITER)
{
    
  ////printf("Start data conversion...\n");
	double *scores;
	short **matrix;
	size_t i,dim,nc,nr,j,e=0;
	scores=NULL;
	nc=asInteger(ncol);
    nr=asInteger(nrow);
  	
    PROTECT(incidence=coerceVector(incidence,INTSXP));
    
    /* Initialize matrix */
    do
        matrix = (short **) R_alloc(nr,sizeof(short*));
    while(matrix==NULL);
    
 	for (i = 0; i<nr;i++)
    {
        do
            matrix[i] = (short *) R_alloc(nc,sizeof(short));
        while(matrix[i]==NULL);
        
    }
  	for (j = 0; j<nc;j++)
		for(i=0;i<nr;i++)
			matrix[i][j] =(short) INTEGER(incidence)[j * nr+i];
	
	for (i=0;i<nr;i++)
	{
		for (j=0;j<nc;j++)
            if(matrix[i][j]==1)e++;
        
	}
    
	//printf("%d %d %d %d \n",nc,nr, asInteger(step),asInteger(max_iter));
	if(asInteger(MAXITER)==0)
	dim=analysis(matrix, nc, nr, &scores, asInteger(step),  asInteger(max_iter),asInteger(verbose));
	else
	dim=analysis_ex(matrix, nc, nr, &scores, asInteger(step),  asInteger(max_iter),asInteger(verbose),asInteger(MAXITER));

 	SEXP res;
    PROTECT(res=allocVector(REALSXP,dim));


  	for (i=0; i<dim; i++)
    {
        REAL(res)[i]=scores[i];

    }
 

   

	SEXP result;
	SEXP names;
	SEXP NN;
	PROTECT(result = allocVector(VECSXP,2));
	PROTECT(names = allocVector(STRSXP, 2));
	PROTECT(NN=allocVector(REALSXP,1));
	REAL(NN)[0]=asInteger(max_iter);
    char *names_c[2] = {"similarity_scores","N"};
    for(i = 0; i < 2; i++)
        SET_STRING_ELT(names,i,mkChar(names_c[i]));
	
	SET_VECTOR_ELT(result, 0, res);
  
	SET_VECTOR_ELT(result, 1, NN );
    
	setAttrib(result, R_NamesSymbol,names);

	UNPROTECT(5);
	return(result);

}




SEXP R_rewire_bipartite (SEXP incidence,SEXP ncol, SEXP nrow,SEXP max_iter, SEXP verbose, SEXP MAXITER)
{
	////printf("Start data conversion...\n");
	short **matrix;
	size_t i,check,nc,nr,j;
	nc=asInteger(ncol);
  	nr=asInteger(nrow);
  	PROTECT(incidence=coerceVector(incidence,INTSXP));
    
    /* Initialize matrix */
    do
        matrix = (short **) R_alloc(nr,sizeof(short*));
    while(matrix==NULL);
    
 	for (i = 0; i<nr;i++)
    {do
        matrix[i] = (short *) R_alloc(nc,sizeof(short));
        while(matrix[i]==NULL);
        
    }
  	for (j = 0; j<nc;j++)
		for(i=0;i<nr;i++)
			matrix[i][j] =(short) INTEGER(incidence)[j * nr+i];
	
    
    
	////printf("Data conversion ok, now I'm computing...\n");
		if(asInteger(MAXITER)!=0)
	check=rewire_bipartite_ex(matrix, nc, nr,  asInteger(max_iter),asInteger(verbose),asInteger(MAXITER));
else
	check=rewire_bipartite(matrix, nc, nr,  asInteger(max_iter),asInteger(verbose));
      if(check==-1  )
		warning("Reached the maximum number admissible of iterations!\n",asInteger(MAXITER)); 
    
    
    
 	SEXP res;
    PROTECT(res=allocVector(INTSXP,nc*nr));
  	for(i=0;i<nc;i++)
		for(j=0;j<nr;j++)
			INTEGER(res)[i*nr+j]=matrix[j][i];

	UNPROTECT(2);
	return(res);
	////printf("Rewiring done!\n");
    
	return 0;
}




SEXP R_rewire_sparse_bipartite (SEXP edges,SEXP ncol, SEXP nrow,SEXP max_iter,SEXP e,SEXP verbose, SEXP MAXITER)
{
	////printf("Start data conversion...\n");
	size_t *from;
	size_t *to;
    
    
	size_t i,check,nc,nr,ne;
	ne=asInteger(e);
	nc=asInteger(ncol);
    nr=asInteger(nrow);
    PROTECT(edges=coerceVector(edges,INTSXP));
    
    /* Initialize matrix */
    do
        from = (size_t *) R_alloc(ne,sizeof(size_t));
    while(from==NULL);
    
    do
        to = (size_t *) R_alloc(ne,sizeof(size_t));
    while(to==NULL);
    
    
    for(i=0;i<ne;i++)
    {
        from[i] =(size_t) INTEGER(edges)[i];
        to[i]		=(size_t) INTEGER(edges)[i+ne];
       // printf("%d %d \n",from[i],to[i]);
        
        
    }
    

	if(asInteger(MAXITER)!=0)
	check=rewire_sparse_bipartite_ex(from,to, nc, nr,  asInteger(max_iter),ne,asInteger(verbose),asInteger(MAXITER));
else
	check=rewire_sparse_bipartite(from,to, nc, nr,  asInteger(max_iter),ne,asInteger(verbose));
      if(check==-1 )
		warning("Reached the maximum number admissible of iterations!\n",asInteger(MAXITER)); 
    
    
    
    
 	SEXP res;
    PROTECT(res=allocVector(INTSXP,ne*2));
    size_t k=0;
  	for(i=0;i<ne;i++)
    {
        INTEGER(res)[k++]=from[i];
        INTEGER(res)[k++]=to[i];
    }

    

	UNPROTECT(2);
	return(res);
	////printf("Rewiring done!\n");
    
	return 0;
}


SEXP R_analysis_undirected (SEXP incidence,SEXP ncol, SEXP nrow, SEXP step, SEXP max_iter,SEXP verbose, SEXP MAXITER)
{

	////printf("Start data conversion...\n");
	double *scores;
	short **matrix;
	size_t i,dim,nc,nr,j,e=0;
	scores=NULL;
	nc=asInteger(ncol);
  nr=asInteger(nrow);
  SEXP incidencetmp;
  	//removed protection, doesn' work with protection
  PROTECT(incidencetmp=coerceVector(incidence,INTSXP));
 //SEXP incidencetmp=coerceVector(incidence,INTSXP);

  /* Initialize matrix */
		do
        matrix = (short **) R_alloc(nr,sizeof(short*));
    while(matrix==NULL);
    
 	for (i = 0; i<nr;i++)
    {
        do
            matrix[i] = (short *) R_alloc(nc,sizeof(short));
        while(matrix[i]==NULL);
        
    }
  	for (j = 0; j<nc;j++)
		for(i=0;i<nr;i++)	    
			matrix[i][j] =(short) INTEGER(incidencetmp)[j * nr+i];
	
	for (i=0;i<nr;i++)
	{	
		for (j=0;j<nc;j++)
				if(matrix[i][j]==1)e++;
			
	}

	////printf("Data conversion ok, now I'm computing...\n");
	if(asInteger(MAXITER)!=0)
	dim=analysis_undirected_ex(matrix, nc, nr, &scores, asInteger(step),  asInteger(max_iter),asInteger(verbose),asInteger(MAXITER));
	else
		dim=analysis_undirected(matrix, nc, nr, &scores, asInteger(step),  asInteger(max_iter),asInteger(verbose));

 
    
 	SEXP res;
    PROTECT(res=allocVector(REALSXP,dim));
  	//restmp=REAL(res);	
	//printf("ok dim= %d ",dim);fflush(stdout);
  	for (i=0; i<dim; i++) 
	  {
			REAL(res)[i]=scores[i];
		//	printf("%lf",(scores)[i] );fflush(stdout);
	  } 	
		

	SEXP result;
	SEXP names;
	SEXP NN;

	PROTECT(result = allocVector(VECSXP, 2));//res,matrix,l,N
	PROTECT(names = allocVector(STRSXP, 2));
	PROTECT(NN=allocVector(REALSXP,1));
	REAL(NN)[0]=asInteger(max_iter);

  char *names_c[2] = {"similarity_scores","N"};
  for(i = 0; i < 2; i++)   
 	 SET_STRING_ELT(names,i,mkChar(names_c[i])); 
	
	SET_VECTOR_ELT(result, 0, res);


	SET_VECTOR_ELT(result, 1, NN );

	setAttrib(result, R_NamesSymbol,names); 

	//UNPROTECT(4);
	UNPROTECT(5);
	return(result);
	////printf("Analysis done!\n");
}




SEXP R_rewire (SEXP incidence,SEXP ncol, SEXP nrow,SEXP max_iter,SEXP verbose, SEXP MAXITER)
{
	////printf("Start data conversion...\n");
	short **matrix;
	size_t i,check,nc,nr,j;
	nc=asInteger(ncol);
  	nr=asInteger(nrow);
  	PROTECT(incidence=coerceVector(incidence,INTSXP));

  /* Initialize matrix */
do
        matrix = (short **) R_alloc(nr,sizeof(short*));
    while(matrix==NULL);
    
 	for (i = 0; i<nr;i++)
    {
        do
            matrix[i] = (short *) R_alloc(nc,sizeof(short));
        while(matrix[i]==NULL);
        
    }
  	for (j = 0; j<nc;j++)
		for(i=0;i<nr;i++)	    
			matrix[i][j] =(short) INTEGER(incidence)[j * nr+i];
	


	////printf("Data conversion ok, now I'm computing...\n");
		if(asInteger(MAXITER)!=0)
	check=rewire_ex(matrix, nc, nr,  asInteger(max_iter),asInteger(verbose),asInteger(MAXITER));
	else
		check=rewire(matrix, nc, nr,  asInteger(max_iter),asInteger(verbose));
      if(check==-1  )
		warning("Reached the maximum number admissible of iterations!\n",asInteger(MAXITER)); 
    



 	SEXP res;
    PROTECT(res=allocVector(INTSXP,nc*nr));
  	for(i=0;i<nc;i++)
		for(j=0;j<nr;j++)
			INTEGER(res)[i*nr+j]=matrix[j][i];

	UNPROTECT(2);
	return(res);
	////printf("Rewiring done!\n");


}


SEXP R_rewire_sparse (SEXP edges,SEXP ncol, SEXP nrow,SEXP max_iter,SEXP e,SEXP verbose, SEXP MAXITER)
{
	////printf("Start data conversion...\n");
	size_t *from;
	size_t *to;
    size_t *degree;

	size_t i,check,nc,nr,ne;
	nc=asInteger(ncol);
  	nr=asInteger(nrow);
    ne=asInteger(e);
  	PROTECT(edges=coerceVector(edges,INTSXP));




do
    from = (size_t *) R_alloc(ne,sizeof(size_t));while(from==NULL);
do
	to = (size_t *) R_alloc(ne,sizeof(size_t));while(to==NULL);
    degree = (size_t *) calloc(nr,sizeof(size_t));
    
    /* Initialize matrix */
    for(i=0;i<ne;i++)
    {
        from[i] =(size_t) INTEGER(edges)[i];
        to[i]		=(size_t) INTEGER(edges)[i+ne];
     //   printf("%d %d \n",from[i],to[i]);
        degree[from[i]]++;
        degree[to[i]]++;

        
    }
	
    
    
	////printf("Data conversion ok, now I'm computing...\n");
if(asInteger(MAXITER)!=0)
			check=rewire_sparse_ex(from,to,degree, nc, nr,  asInteger(max_iter),ne,asInteger(verbose),asInteger(MAXITER));
else
check=rewire_sparse(from,to,degree, nc, nr,  asInteger(max_iter),ne,asInteger(verbose));
      if(check==-1  )
		warning("Reached the maximum number admissible of iterations!\n",asInteger(MAXITER)); 
    
    
    SEXP res;
    PROTECT(res=allocVector(INTSXP,ne*2));
    size_t k=0;
  	for(i=0;i<ne;i++)
    {
        INTEGER(res)[k++]=from[i];
        INTEGER(res)[k++]=to[i];
    }

	UNPROTECT(2);
	return(res);
	////printf("Rewiring done!\n");
    
	return 0;
}



