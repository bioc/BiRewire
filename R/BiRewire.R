# This code is written by Andrea Gobbi  <gobbi.andrea@mail.com>
# 2015
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
# You should have received a copy of the GNU General Public licenses
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
birewire.analysis<- function(incidence, step=10, max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
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


		if(exact==TRUE)
    	{
					if( max.iter=="n")
						max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
					MAXITER_MUL=MAXITER_MUL*max.iter
  				result<-.Call("R_analysis", incidence,nc,nr,as.numeric(step),as.numeric(max.iter),verbose,MAXITER_MUL+1)
					result$N=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
			}

		else
			{
					if( max.iter=="n")
    				max.iter=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy) )  
  				result<-.Call("R_analysis", incidence,nc,nr,as.numeric(step),as.numeric(max.iter),verbose,0)
    				result$N=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy)   )


			}
  #C function is called!
	#cat(c(nc,nr,step,max.iter,"\n"))



  return( result)
}
##Performs the rewiring algorithm max.iter times
birewire.rewire.bipartite<- function(incidence,  max.iter="n", accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
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

	
		if(exact==TRUE)
    	{
					if( max.iter=="n")
						max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
					MAXITER_MUL=MAXITER_MUL*max.iter
  				result<-.Call("R_rewire_bipartite", incidence,nc , nr, as.numeric( max.iter),verbose,MAXITER_MUL+1)

			}

		else
			{
					if( max.iter=="n")
    				max.iter=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy) )  
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
birewire.rewire.bipartite.and.projections<-function(graph,step=10,max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10)
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
    max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
  sc=c()
  sc1=c()
  sc2=c()
  sc[1]=sc1[1]=sc2[1]=1
  mm=m
  for(i in 1:(max.iter/step))
  {
    mm=birewire.rewire.bipartite(incidence=mm,max.iter=step,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=FALSE)
    p=birewire.bipartite.from.incidence(directed=FALSE,matrix=mm,reverse=FALSE)
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
birewire.rewire.sparse.bipartite<- function(graph,  max.iter="n", accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
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
					  max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
					MAXITER_MUL=MAXITER_MUL*max.iter
 					 result<-.Call("R_rewire_sparse_bipartite", edges,nc , nr, as.numeric( max.iter),e,verbose,MAXITER_MUL+1)
  

			}

		else
			{
					if( max.iter=="n")
    				max.iter=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy) )  
 					 result<-.Call("R_rewire_sparse_bipartite", edges,nc , nr, as.numeric( max.iter),e,verbose,0)



			}
  
  
  gg<-graph.bipartite(edges=result,types=V(g)$type,directed=TRUE)
  
  
  return(gg)
}

birewire.analysis.undirected<- function(adjacency, step=10, max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE) 
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



		if(exact==TRUE)
    	{
						if( max.iter=="n")
   						 max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
						MAXITER_MUL=MAXITER_MUL*max.iter
						result<-.Call("R_analysis_undirected", adjacency,n,n,as.numeric(step),as.numeric(max.iter),as.numeric(verbose),MAXITER_MUL)
						result$N=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
  

			}

		else
			{

					if(max.iter=="n")
			 			max.iter=(e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy)
					result<-.Call("R_analysis_undirected", adjacency,n,n,as.numeric(step),as.numeric(max.iter),as.numeric(verbose),0)
					result$N=(e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy)


			}


		return( result)
	}



birewire.rewire<- function(adjacency,  max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
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

	if(exact==TRUE)
    	{
						if( max.iter=="n" )
    					max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
						MAXITER_MUL=MAXITER_MUL*max.iter
 						result<-.Call("R_rewire", adjacency,n , n, as.numeric( max.iter),verbose,MAXITER_MUL+1)


			}

		else
			{

					if(max.iter=="n")
			 			max.iter=(e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy)
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


birewire.rewire.sparse<- function(graph,  max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
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
        
        

	if(exact==TRUE)
    	{
						if( max.iter=="n")
  				 	 max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
						MAXITER_MUL=MAXITER_MUL*max.iter
 						result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,MAXITER_MUL+1)


			}

		else
			{

        		if(max.iter=="n")
			 		max.iter=(e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy)
				
 				result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,0)

			}

		
        gg<-graph(edges=result+1,directed=FALSE,n=n)
		return( gg)
	}











#######V
##Given a dsg (a list of two incidnece matrix), this routine sample randomly K network preserving the degrees. The inputs are the same of birewire.rewire.bipartite
birewire.sampler.dsg<-function(dsg,K,path,delimitators=list(negative='-',positive='+'),exact=FALSE,verbose=TRUE, max.iter.pos='n',max.iter.neg='n', accuracy=0.00005,MAXITER_MUL=10)
	{

  
		if(!is.list(dsg) | is.null(dsg[['positive']]) | is.null(dsg[['negative']]) )
			    {
			    	stop("The input must be a dsg object (see References) \n")
			    	return(0)
			    }

		if(!file.exists(path))
  					{
    					dir.create(path) 
 				 	}

		n=ceiling(K/1000)
		NNET=1000
		if(n==0)
			n=1
		for( i in 1:n)
			{
				if(K-1000*i<0)
					NNET=K-1000*(i-1)			
  
			  	print(paste('Filling directory n.',i,'with',NNET,'randomised versions of the given dsg.'))
    			PATH<-paste(path,'/',i,'/',sep='')
				if(!file.exists(PATH))
					{
      					dir.create(PATH)
    				}
    				for(j in 1:NNET)
    					{
    						dsg=birewire.rewire.dsg(dsg=dsg,delimitators=delimitators,exact=exact,path=paste(PATH,'network_',(i-1)*1000+j,'.sif',sep=''),
    							verbose=verbose,max.iter.pos=max.iter.pos,max.iter.neg=max.iter.neg, accuracy=accuracy,MAXITER_MUL=MAXITER_MUL)

    					}	

			}





	}


##sililat to above function
birewire.sampler.bipartite<-function(incidence,K,path,max.iter="n", accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE,write.sparse=TRUE)
{


		if(!file.exists(path))
  					{
    					dir.create(path) 
 				 	}

		n=ceiling(K/300)
		##NB 300 perche' non voglio piu' di 1000 file per cartella
		NFILES=300
		if(write.sparse==F)
			NFILES=1000
		NNET=NFILES

		for( i in 1:n)
			{
					if(K-NFILES*i<0)
					NNET=K-NFILES*(i-1)
  
			  	print(paste('Filling directory n.',i,'with',NNET,'randomised versions of the given dsg.'))
    			PATH<-paste(path,'/',i,'/',sep='')
				if(!file.exists(PATH))
					{
      					dir.create(PATH)
    				}
    				for(j in 1:NNET)
    					{
    						incidence=birewire.rewire.bipartite(incidence=incidence,  max.iter=max.iter, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,
    							exact=exact)
    						if(write.sparse)
    							{
    									write_stm_CLUTO(as.simple_sparse_array(as.matrix(incidence)),file=paste(PATH,'network_',(i-1)*1000+j,sep=''))
    							}else
    							{
										write.table(incidence,file=paste(PATH,'network_',(i-1)*1000+j,sep=''),append=F)

    							}

    					}	

			}




	

}

birewire.rewire.dsg<-function(dsg,exact=FALSE,verbose=1,max.iter.pos='n',max.iter.neg='n',accuracy=0.00005,MAXITER_MUL=10,path=NULL,delimitators=list(positive='+',negative= '-'))
{


		if(!is.list(dsg) | is.null(dsg[['positive']]) | is.null(dsg[['negative']]) )
			    {
			    	stop("The input must be a dsg object (see References) \n")
			    	return(0)
			    }


	incidence_pos=dsg[["positive"]]
	incidence_neg=dsg[["negative"]]
	incidence_pos=birewire.rewire.bipartite(incidence=incidence_pos,  max.iter=max.iter.pos, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
    incidence_neg=birewire.rewire.bipartite(incidence=incidence_neg,  max.iter=max.iter.neg, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
	dsg=list(positive=incidence_pos,negative=incidence_neg)	
	if(!is.null(path))
		{

			birewire.save.dsg(g=birewire.build.dsg(dsg,delimitators),file=path)
				
		}
	return(dsg)
}


##INTERNAL ROUTINES
simplify.table<-function(df)
	{
		return(df[which(rowSums(df)>0),which(colSums(df)>0)])

	}
##INTERNAL ROUTINES
get.data.frame.from.incidence<-function(table,sign)
{


source=rownames(table)
target=colnames(table)
index=which(table>0,arr.ind=TRUE)
df=data.frame(source=source[index[,1]],sign=sign,target=target[index[,2]])
return(df)

}

##from a sif file, the routine generates the negative and positive incidence matrix
birewire.induced.bipartite<-function(g,delimitators=list(negative='-',positive='+'))

{
	if(is.null(delimitators[['positive']]) | is.null(delimitators[['negative']]) )
			    {
			    	stop("Problem with delimitators (see References) \n")
			    	return(0)
			    }


	g=as.data.frame(g)
	dsg=list()
	g_p=g[g[,2]==delimitators[['positive']],c(1,3)]
	g_n=g[g[,2]==delimitators[['negative']],c(1,3)]

	positive=as.data.frame.matrix(table(g_p))

	negative=as.data.frame.matrix(table(g_n))
	dsg[['positive']]=simplify.table(positive)
	
	dsg[['negative']]=simplify.table(negative)
	
return(dsg)

}	


##inverse of the function above
birewire.build.dsg<-function(dsg,delimitators=list(negative='-',positive='+'))
{


	if(!is.list(dsg) | is.null(dsg[['positive']]) | is.null(dsg[['negative']]) )
			    {
			    	stop("The input must be a dsg object (see References) \n")
			    	return(0)
			    }
	positive=dsg[['positive']]
	negative=dsg[['negative']]

	g_p=get.data.frame.from.incidence(positive,delimitators[['positive']])
	g_n=get.data.frame.from.incidence(negative,delimitators[['negative']])
	g=rbind(g_p,g_n)
	return(g)
}

birewire.load.dsg<-function(path)
	{
		

		return(unique(read.table(path)))


	}
birewire.save.dsg<-function(g,file)
	{

		
		write.table(g,file,col.names=FALSE,row.names=FALSE,quote=FALSE)

		


	}

##jaccard index for dsg
	birewire.similarity.dsg<-function(m1,m2)
{
	
 x=sum(m1[['positive']]*m2[['positive']]) +sum(m1[['negative']]*m2[['negative']] )
e=sum(m1[['positive']])+sum(m2[['positive']])+sum(m1[['negative']])+sum(m2[['negative']])
  return( x/(e-x))
  
}

