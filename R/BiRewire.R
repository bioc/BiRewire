# This code is written by Andrea Gobbi  <gobbi.andrea@mail.com>
# 2013
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 



# Create a bipartite igraph graph from its incidence matrix
# it could be directed and the two classes can be reversed
birewire.bipartite.from.incidence<-function(matrix,directed=FALSE,reverse=FALSE)
{ 

	edges=which(matrix==1,arr.ind=TRUE)

		edges[,2]=	edges[,2]+nrow(matrix)
  cl1=rep(T,nrow(matrix))
  cl2=rep(F,ncol(matrix))
  if(reverse)
  {
    cl1=!cl1
    cl2=!cl2
  }
  types=c(cl1,cl2)
g=graph.bipartite(types,as.vector(t(edges)),directed)
if(length(rownames(matrix))*length(colnames(matrix))!=0)
V(g)$label=c(rownames(matrix),colnames(matrix))
return(g)
}
##Pertforms the analysis of a bipartite graph (see the manual for more details)
##incidence= the incidence matrix of the graph
##step=fraquancy for which the score is calculated (if "n" a bound is computed)
##max.iter=maximum number of successfull rewiring steps
##accuracy=level of accuracy
##verbose= print execution bar during the process 
##MAXITER_MUL= MAXITER_MUL* max.iter indicates the maximum number of real iteration
birewire.analysis<- function(incidence, step=10, max.iter="n",accuracy=1,verbose=TRUE,MAXITER_MUL=10,exact=F)
{
  
if(!is.matrix(incidence) && !is.data.frame(incidence))
	{
    stop("The input must be a data.frame or a matrix object \n")
    return (0)
  }

		if(is.data.frame(incidence))
		{	
			incidence=as.matrix(incidence)
		}
	if(verbose)
    verbose=1
  else
    verbose=0
  nc=as.numeric(ncol(incidence))
  nr=as.numeric(nrow(incidence))
  t=nc*nr
	

	

	if(is.integer(incidence)|| is.logical(incidence))
		{
			incidence=matrix(as.double(incidence),ncol=nc)
		}


  e=sum(incidence)


		if(exact==T)
    	{
					if( max.iter=="n")
						max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
					MAXITER_MUL=MAXITER_MUL*max.iter
  				result<-.Call("R_analysis", incidence,nc,nr,as.numeric(step),as.numeric(max.iter),verbose,MAXITER_MUL+1)
					result$N=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
			}

		else
			{
					if( max.iter=="n")
    				max.iter=floor((e/(2-2*e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )  )
  				result<-.Call("R_analysis", incidence,nc,nr,as.numeric(step),as.numeric(max.iter),verbose,0)
    				result$N=floor((e/(2-2*e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )  )


			}
  #C function is called!
	#cat(c(nc,nr,step,max.iter,"\n"))



  return( result)
}
##Performs the rewiring algorithm max.iter times
birewire.rewire.bipartite<- function(incidence,  max.iter="n", accuracy=1,verbose=TRUE,MAXITER_MUL=10,exact=F)
{
  

if(!is.matrix(incidence) && !is.data.frame(incidence)  && !is.igraph(incidence))
	{
    stop("The input must be a data.frame or a matrix object or an igraph bipartite graph \n")
    return (0)

  }
if(is.igraph(incidence))
	return(birewire.rewire.sparse.bipartite(incidence,  max.iter, accuracy,verbose,MAXITER_MUL,exact))

	dataframe=FALSE
		if(is.data.frame(incidence))
		{	dataframe=TRUE
			incidence=as.matrix(incidence)
		}

  if(verbose)
    verbose=1
  else
    verbose=0
  nc=as.numeric(ncol(incidence))
  nr=as.numeric(nrow(incidence))
  t=nc*nr
	rnames=rownames(incidence)
	cnames=colnames(incidence)
	out="double"
	if(is.integer(incidence))
		out="integer"
	if(is.logical(incidence))
		out="logical"

	if(is.integer(incidence)|| is.logical(incidence))
		{
			incidence=matrix(as.double(incidence),ncol=nc)
		}


  e=sum(incidence)

	
		if(exact==T)
    	{
					if( max.iter=="n")
						max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
					MAXITER_MUL=MAXITER_MUL*max.iter
  				result<-.Call("R_rewire_bipartite", incidence,nc , nr, as.numeric( max.iter),verbose,MAXITER_MUL+1)

			}

		else
			{
					if( max.iter=="n")
    				max.iter=floor((e/(2-2*e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )  )
 					result<-.Call("R_rewire_bipartite", incidence,nc , nr, as.numeric( max.iter),verbose,0)



			}

 	if(out=="double")
 		result=matrix(result,ncol=ncol(incidence))
	else
 			result=matrix(as.integer(result),ncol=ncol(incidence))
	if(out=="logical")
 				result=matrix(as.logical(result),ncol=ncol(incidence))
		if(dataframe)
			result=as.data.frame(result)
	colnames(result)=cnames
	rownames(result)=rnames
  return( result)
}

##SA algorithm for monitoring the behaviour of the two natural projections (Slow)
birewire.rewire.bipartite.and.projections<-function(graph,step=10,max.iter="n",accuracy=1,verbose=TRUE,MAXITER_MUL=10)
{

  g=graph
  if(is.bipartite(g)==FALSE)
  {
    stop("Not a bipartite graph \n")
    return (0)
  }
  g1= bipartite.projection(g,multiplicity=FALSE)$proj1
  m1=get.adjacency(graph=g1,sparse=FALSE)
  g2= bipartite.projection(g,multiplicity=FALSE)$proj2
  m2=get.adjacency(graph=g2,sparse=FALSE)
  m=get.incidence(g)
  nr=nrow(m)
  nc=ncol(m)



  e=sum(m)
  t=nr*nc
  if( max.iter=="n")
    max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
  sc=c()
  sc1=c()
  sc2=c()
  sc[1]=sc1[1]=sc2[1]=1
  mm=m
  for(i in 1:(max.iter/step))
  {
    mm=birewire.rewire.bipartite(incidence=mm,max.iter=step,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=T)
    p=birewire.bipartite.from.incidence(directed=F,matrix=mm,reverse=FALSE)
    sc[i+1]=birewire.similarity(m,mm)
    p=bipartite.projection(graph=p,multiplicity=FALSE)
    tmp1=get.adjacency(graph=p$proj1,sparse=FALSE)
    tmp2=get.adjacency(graph=p$proj2,sparse=FALSE)
    if(dim(tmp1)[1]==dim(m1)[1])
    {
      sc1[i+1]=birewire.similarity(m1,  tmp1 )
      sc2[i+1]=birewire.similarity(m2,  tmp2)
    }
    else
    {
      sc1[i+1]=birewire.similarity(m2,  tmp1)
      sc2[i+1]=birewire.similarity(m1,  tmp2)
    }
  }
  result=list()
  result$proj1=g1
  result$proj2=g2
  result$similarity_scores.proj1=sc1
  result$similarity_scores.proj2=sc2
  result$similarity_scores=sc
  result$rewired.proj1=p$proj1
  result$rewired.proj2=p$proj2
  result$rewired=birewire.bipartite.from.incidence(directed=F,matrix=mm,reverse=FALSE)
  return( result )
}
#similarity score (jaccard and hamming)
birewire.similarity<-function(m1,m2)
{
  if(dim(m2)[1]!=dim(m1)[1])
    m1=t(m1)
  return( sum( m1*m2)/sum(m1+m2-m1*m2))
  
}
birewire.rewire.sparse.bipartite<- function(graph,  max.iter="n", accuracy=1,verbose=TRUE,MAXITER_MUL=10,exact=F)
{
  if(verbose)
    verbose=1
  else
    verbose=0
  g=graph
  if(is.bipartite(g)==FALSE)
  {
    stop("Not a bipartite graph \n")
    return (0)
  }

  
  
  e=length(E(g))
  edges= get.edgelist(names=FALSE,g)
  edges=edges[order(edges[,1]),]
  
  nr=length(unique(edges[,1]))
  nc=length(V(g))-nr
  t=nc*nr
		if(exact==T)
    	{
					if( max.iter=="n")
					  max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
					MAXITER_MUL=MAXITER_MUL*max.iter
 					 result<-.Call("R_rewire_sparse_bipartite", edges,nc , nr, as.numeric( max.iter),e,verbose,MAXITER_MUL+1)
  

			}

		else
			{
					if( max.iter=="n")
    				max.iter=floor((e/(2-2*e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )  )
 					 result<-.Call("R_rewire_sparse_bipartite", edges,nc , nr, as.numeric( max.iter),e,verbose,0)



			}
  
  
  gg<-graph.bipartite(edges=result,types=V(g)$type,directed=TRUE)
  
  
  return(gg)
}

birewire.analysis.undirected<- function(adjacency, step=10, max.iter="n",accuracy=1,verbose=TRUE,MAXITER_MUL=10,exact=F) 
	{

if(!is.matrix(adjacency) && !is.data.frame(adjacency))
	{
    stop("The input must be a data.frame or a matrix object \n")
    return (0)
  }

		if(is.data.frame(adjacency))
			adjacency=as.matrix(adjacency)
				
			
		 if(verbose)
    verbose=1
  else
    verbose=0
		n=as.numeric(ncol(adjacency))
		if(is.integer(adjacency)|| is.logical(adjacency))
		{
			adjacency=matrix(as.double(adjacency),ncol=n)
		}
		t=n^2/2
		e=sum(adjacency)/2
		d=e/t



		if(exact==T)
    	{
						if( max.iter=="n")
   						 max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
						MAXITER_MUL=MAXITER_MUL*max.iter
						result<-.Call("R_analysis_undirected", adjacency,n,n,as.numeric(step),as.numeric(max.iter),as.numeric(verbose),MAXITER_MUL)
						result$N=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
  

			}

		else
			{

					x=(accuracy*n^4 )/(n^4*e - 2*n^2*e^2)		
					base=((e - 2)*n^6 - 4*n^4*e + 24*n^2*e^2 - 16*e^3)*e^(-1)/n^6
					if( max.iter=="n" & accuracy!=1)
							 max.iter=floor(log(x=x,base= base))	
					if( max.iter=="n" & accuracy==1)
			 				max.iter=(e/(2*d^3-6*d^2+2*d+2))*log(e*(1-d))
					result<-.Call("R_analysis_undirected", adjacency,n,n,as.numeric(step),as.numeric(max.iter),as.numeric(verbose),0)
					result$N=floor(log(x=x,base= base))
					if( max.iter=="n" & accuracy==1)
						 result$N=(e/(2*d^3-6*d^2+2*d+2))*log(e*(1-d))


			}


		return( result)
	}



birewire.rewire<- function(adjacency,  max.iter="n",accuracy=1,verbose=TRUE,MAXITER_MUL=10,exact=F)
	{ 

if(!is.matrix(adjacency) && !is.data.frame(adjacency) && !is.igraph(adjacency) )
	{
    stop("The input must be a data.frame or a matrix object \n")
    return (0)
  }
		dataframe=FALSE
		if(is.data.frame(adjacency))
			{
				adjacency=as.matrix(adjacency)
				dataframe=TRUE
			}
		if(is.igraph(adjacency))
			return(birewire.rewire.sparse(adjacency,  max.iter,accuracy,verbose,MAXITER_MUL,exact))

	if(verbose)
    verbose=1
  else
    verbose=0
 		n=as.numeric(ncol(adjacency))

		rnames=rownames(adjacency)
		cnames=colnames(adjacency)
		out="double"
		if(is.integer(adjacency))
			out="integer"
		if(is.logical(adjacency))
			out="logical"

		if(is.integer(adjacency)|| is.logical(adjacency))
		{
			adjacency=matrix(as.double(adjacency),ncol=n)
		}


		t=n^2/2
		e=sum(adjacency)/2
		d=e/t

	if(exact==T)
    	{
						if( max.iter=="n" )
    					max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
						MAXITER_MUL=MAXITER_MUL*max.iter
 						result<-.Call("R_rewire", adjacency,n , n, as.numeric( max.iter),verbose,MAXITER_MUL+1)


			}

		else
			{

					x=(accuracy*n^4)/(n^4*e  - 2*n^2*e^2)
					base=((e - 2)*n^6 - 4*n^4*e + 24*n^2*e^2 - 16*e^3)*e^(-1)/n^6
					if( max.iter=="n" & accuracy!=1)
			 			max.iter=floor(log(x=x,base= base))	
					if( max.iter=="n" & accuracy==1)
						 max.iter=(e/(2*d^3-6*d^2+2*d+2))*log(e*(1-d))
 					result<-.Call("R_rewire", adjacency,n , n, as.numeric( max.iter),verbose,0)


			}


 	if(out=="double")
 		result=matrix(result,ncol=ncol(adjacency))
	else
 			result=matrix(as.integer(result),ncol=ncol(adjacency))
	if(out=="logical")
 				result=matrix(as.logical(result),ncol=ncol(adjacency))
		if(dataframe)
			result=as.data.frame(result)
		colnames(result)=cnames
		rownames(result)=rnames

		return( result)
	}


birewire.rewire.sparse<- function(graph,  max.iter="n",accuracy=1,verbose=TRUE,MAXITER_MUL=10,exact=F)
	{ if(verbose)
    verbose=1
  else
    verbose=0
 		n=as.numeric(length(V(graph)))
        
		t=n^2/2
		e=as.numeric(length(E(graph)))

		d=e/t
        edges= get.edgelist(graph)
        edges[ , c(1,2)] <- edges[ , c(2,1)]
        edges=edges[order(edges[,1]),]-1
        
        

	if(exact==T)
    	{
						if( max.iter=="n")
  				 	 max.iter=floor((e*(1-e/t)) *log(x= e/accuracy-e*e/(accuracy*t) )/2  )
						MAXITER_MUL=MAXITER_MUL*max.iter
 						result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,MAXITER_MUL+1)


			}

		else
			{

							x=(accuracy*n^4)/(n^4*e  - 2*n^2*e^2)
        

		base=((e - 2)*n^6 - 4*n^4*e + 24*n^2*e^2 - 16*e^3)*e^(-1)/n^6
		if( max.iter=="n" & accuracy!=1)
			 max.iter=floor(log(x=x,base= base))	
		if( max.iter=="n" & accuracy==1)
			 max.iter=(e/(2*d^3-6*d^2+2*d+2))*log(e*(1-d))
 		result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,0)


			}

		
        gg<-graph(edges=result+1,directed=FALSE,n=n)
		return( gg)
	}

