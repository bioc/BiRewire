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
# it could be directed and the two classes
birewire.bipartite.from.incidence<-function(matrix,directed=FALSE)
{ 
	return(graph.incidence(matrix,directed))
}
##Pertforms the analysis of a bipartite graph (see the manual for more details)
##incidence= the incidence matrix of the graph
##step=fraquancy for which the score is calculated (if "n" a bound is computed)
##max.iter=maximum number of successfull rewiring steps
##accuracy=level of accuracy
##verbose= print execution bar during the process 
##MAXITER_MUL= MAXITER_MUL* max.iter indicates the maximum number of real iteration
birewire.analysis.bipartite<- function(incidence, step=10, max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE,n.networks=50,display=TRUE)
{
  
	if(is.null(incidence))
		return( list(N=0,data=NULL))

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
	RES=NULL
	if(exact==TRUE)
    	{
			if( max.iter=="n")
				max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
			MAXITER_MUL=MAXITER_MUL*as.numeric(max.iter)
		}else
		{
			if( max.iter=="n")
				max.iter=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy) )  
		}

	for(i in 1:n.networks)
	{
		if(exact)
			{
			result<-.Call("R_analysis", incidence,nc,nr,as.numeric(step),as.numeric(max.iter),verbose,MAXITER_MUL+1)
			result$N=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
		}
		else
		{
				result<-.Call("R_analysis", incidence,nc,nr,as.numeric(step),as.numeric(max.iter),verbose,0)
			result$N=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy)   )


		}
	RES=rbind(RES,result$similarity_scores)

	}
	if(display)
	{
		mean=colMeans(RES)
		std=apply(RES,2,sd)
		
		sup=mean+qt(.975,nrow(RES)-1)*std/sqrt(nrow(RES))
		inf=mean-qt(.975,nrow(RES)-1)*std/sqrt(nrow(RES))
		par(mfrow=c(2,1))
		x=seq(1,length.out=length(mean))
		plot(step*x,mean,type= 'n',col='blue',lwd=2,main="Jaccard index (JI) over time",xlab="Switching steps",ylab='Jaccard Index')
		polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey80', border = NA)
		lines(step*x,mean,col='blue',lwd=2)
		abline(v=result$N,col= 'red')
		legend("topright",
				col=c('blue','grey80','red'),
				lwd=c(2,10,1),
				legend=c( "Mean JI","C.I.","Bound")
				)
		plot(step*x,mean,type= 'n',col='blue',lwd=2,main="Jaccard index (JI) over time (log-log scale)",log='xy',xlab="Switching steps",ylab='Jaccard Index')
		polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey80', border = NA)
		lines(step*x,mean,col='blue',lwd=2)
		abline(v=result$N,col= 'red')
		legend("bottomleft",
				col=c('blue','grey80','red'),
				lwd=c(2,10,1),
				legend=c( "Mean JI","C.I.","Bound")
				)
	}
	return( list(N=result$N,data=RES))
}
##Performs the rewiring algorithm max.iter times
birewire.rewire.bipartite<- function(incidence,  max.iter="n", accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
{
  
if(is.null(incidence))
	return(NULL)
if(!is.matrix(incidence) && !is.data.frame(incidence)  && !is.igraph(incidence))
	{
    stop("The input must be a data.frame or a matrix object or an igraph bipartite graph \n")
    return (0)

  }
if(is.igraph(incidence))
	return(birewire.rewire.sparse.bipartite(incidence,  max.iter= max.iter, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact))	
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
					MAXITER_MUL=MAXITER_MUL*as.numeric(max.iter)
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
    p=birewire.bipartite.from.incidence(directed=FALSE,matrix=mm)
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
  result$rewired=birewire.bipartite.from.incidence(directed=F,matrix=mm)
  return( result )
}
#similarity jaccard score 
birewire.similarity<-function(m1,m2)
{
	#print("mi rompo")
	if(is.igraph(m1))
	{
		#if(is.bipartite(m1))
		#{
		#	m1=get.incidence(m1,sparse=F)
		#	m2=get.incidence(m2,sparse=F)
		#}else
		#{
		#	m1=get.adjacency(m1,sparse=F)
		#m2=get.adjacency(m2,sparse=F)
		#}
		e=length(E(m1))
		x=length(E(graph.intersection(as.undirected(m1),as.undirected(m2))))
		return(x/(2*e-x))

	}else{
	
	if(dim(m2)[1]!=dim(m1)[1])
		m1=t(m1)
	return( sum( m1*m2)/sum(m1+m2-m1*m2))
	}
  
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
  ##NB sholud have all id from 1 to nnodes
  edges=get.edgelist(names=FALSE,g)

  edges=edges[order(edges[,1]),]-1
  
  nr=length(unique(edges[,1]))
  nc=length(V(g))-nr
  t=nc*nr
  names=V(g)$name
  types=V(g)$type
		if(exact==TRUE)
    	{
					if( max.iter=="n")
					  max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
					MAXITER_MUL=MAXITER_MUL*as.numeric(max.iter)
 					 result<-.Call("R_rewire_sparse_bipartite", edges,nc , nr, as.numeric( max.iter),e,verbose,MAXITER_MUL+1)
  

			}

		else
			{
					if( max.iter=="n")
    				max.iter=ceiling((e/(2-2*e/t)) *log(x=(1-e/t)/accuracy) )  

 					 result<-.Call("R_rewire_sparse_bipartite", edges,nc , nr, as.numeric( max.iter),e,verbose,0)



			}
	#print(edges+1)
  	#print(result+1)
	if(!is.null(names))
	{
		gg<-graph.edgelist(t(matrix(names[result+1],nrow=2)))
		#gg= graph.bipartite( types,names[result+1])
	}else
	{
		gg<-graph.edgelist(t(matrix(result+1,nrow=2)))
		#gg= graph.bipartite( types,result+1)

	}
	V(gg)$type=types
  	
  return(gg)
}

birewire.analysis.undirected<- function(adjacency, step=10, max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE,n.networks=50,display=TRUE) 
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
	RES=NULL
	if(exact==TRUE)
		{
			if( max.iter=="n")
			max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
			MAXITER_MUL=MAXITER_MUL*as.numeric(max.iter)
			}else
			{
				if(max.iter=="n")
			max.iter=ceiling((e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy))
	
			}
	for( i in 1:n.networks)
	{
		if(exact==TRUE)
		{
			result<-.Call("R_analysis_undirected", adjacency,n,n,as.numeric(step),as.numeric(max.iter),as.numeric(verbose),MAXITER_MUL)
			result$N=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
		}

		else
		{
			result<-.Call("R_analysis_undirected", adjacency,n,n,as.numeric(step),as.numeric(max.iter),as.numeric(verbose),0)
			result$N=(e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy)
		}
		RES=rbind(RES,result$similarity_scores)

	}
	if(display)
	{
		mean=colMeans(RES)
		std=apply(RES,2,sd)
		sup=mean+qt(.975,nrow(RES)-1)*std/sqrt(nrow(RES))
		inf=mean-qt(.975,nrow(RES)-1)*std/sqrt(nrow(RES))
		par(mfrow=c(2,1))
		x=seq(1,length.out=length(mean))
		plot(step*x,mean,type= 'n',col='blue',lwd=2,main="Jaccard index (JI) over time",xlab="Switching steps",ylab='Jaccard Index')
		polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey80', border = NA)
		lines(step*x,mean,col='blue',lwd=2)
		abline(v=result$N,col= 'red')
		legend("topright",
				col=c('blue','grey80','red'),
				lwd=c(2,10,1),
				legend=c( "Mean JI","C.I.","Bound")
				)
		plot(step*x,mean,type= 'n',col='blue',lwd=2,main="Jaccard index (JI) over time (log-log scale)",log='xy',xlab="Switching steps",ylab='Jaccard Index')
		polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey80', border = NA)
		lines(step*x,mean,col='blue',lwd=2)
		abline(v=result$N,col= 'red')
		legend("bottomleft",
				col=c('blue','grey80','red'),
				lwd=c(2,10,1),
				legend=c( "Mean JI","C.I.","Bound")
				)

	}
	return( list(N=result$N,data=RES))
}



birewire.rewire.undirected<- function(adjacency,  max.iter="n",accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
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
						MAXITER_MUL=MAXITER_MUL*as.numeric(max.iter)
 						result<-.Call("R_rewire", adjacency,n , n, as.numeric( max.iter),verbose,MAXITER_MUL+1)


			}

		else
			{

					if(max.iter=="n")
			 			max.iter=ceiling((e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy))
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
	{ 
	if(verbose)
		verbose=1
	else
		verbose=0
	n=as.numeric(length(V(graph)))
	t=n^2/2
	e=as.numeric(length(E(graph)))
	d=e/t
	edges= get.edgelist(graph,names=F)
	edges=edges[order(edges[,1]),]-1
	names=V(graph)$name


	if(exact==TRUE)
	{
	if( max.iter=="n")
		 max.iter=ceiling((e*(1-e/t)) *log(x=(1-e/t)/accuracy) /2  )
	MAXITER_MUL=MAXITER_MUL*as.numeric(max.iter)
	result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,MAXITER_MUL+1)


	}

	else
	{

	if(max.iter=="n")
	max.iter=ceiling((e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy))
	result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,0)

	}

	if(!is.null(names))
	{
		gg<-graph.edgelist(t(matrix(names[result+1],nrow=2)),directed=FALSE)
	}else
	{
				gg<-graph.edgelist(t(matrix(result+1,nrow=2)),directed=FALSE)

	}
	return( gg)
	}







birewire.sampler.bipartite<-function(incidence,K,path,max.iter="n", accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE,write.sparse=TRUE)
{
	if(K>=100000)
		{

			stop("I can not create more than 1000 subfolders but if you need it contact the manteiner.")
		}

		if(!file.exists(path))
  					{
    					dir.create(path) 
 				 	}

		
		##NB 300 perche' non voglio piu' di 1000 file per cartella
		NFILES=300
		if(write.sparse==F)
			NFILES=1000
		NNET=NFILES
		n=ceiling(K/NFILES)
		for( i in 1:n)
			{
					if(K-NFILES*i<0)
					NNET=K-NFILES*(i-1)
  
			  	print(paste('Filling directory n.',i,'with',NNET,'randomised versions of the given bipartite.'))
    			PATH<-paste(path,'/',i,'/',sep='')
				if(!file.exists(PATH))
					{
      					dir.create(PATH)
    				}
    				for(j in 1:NNET)
    					{
    						incidence=birewire.rewire.bipartite(incidence=incidence,  max.iter=max.iter, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
    						if(is.igraph(incidence))
								{
									if(write.sparse)
									{
										write_stm_CLUTO(as.simple_sparse_array(as.matrix(get.incidence(incidence,names=TRUE,sparse=FALSE))),file=paste(PATH,'network_',(i-1)*1000+j,sep=''))
									}else
									{
										write.table(get.incidence(incidence,names=TRUE,sparse=FALSE),file=paste(PATH,'network_',(i-1)*1000+j,sep=''),append=F)
									}
								}else
								{
									if(write.sparse)
									{
										write_stm_CLUTO(as.simple_sparse_array(as.matrix(incidence,names=TRUE,names=TRUE)),file=paste(PATH,'network_',(i-1)*1000+j,sep=''))
									}else
									{
										write.table(incidence,file=paste(PATH,'network_',(i-1)*1000+j,sep=''),append=FALSE)
									}
								}

    					}	

			}




	

}



birewire.sampler.undirected<-function(adjacency,K,path,max.iter="n", accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE,write.sparse=TRUE)
{
	if(K>=100000)
		{

			stop("I can not create more than 1000 subfolders but if you need it contact the manteiner.")
		}

		if(!file.exists(path))
  					{
    					dir.create(path) 
 				 	}

		##NB 300 perche' non voglio piu' di 1000 file per cartella
		NFILES=300
		if(write.sparse==F)
			NFILES=1000
		NNET=NFILES
		n=ceiling(K/NFILES)

		for( i in 1:n)
			{
					if(K-NFILES*i<0)
					NNET=K-NFILES*(i-1)
  
			  	print(paste('Filling directory n.',i,'with',NNET,'randomised versions of the given undirected graph.'))
    			PATH<-paste(path,'/',i,'/',sep='')
				if(!file.exists(PATH))
					{
      					dir.create(PATH)
    				}
    				for(j in 1:NNET)
    					{
    						adjacency=birewire.rewire.undirected(adjacency=adjacency,  max.iter=max.iter, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
    						if(is.igraph(adjacency))
								{
									if(write.sparse)
									{
										write_stm_CLUTO(as.simple_sparse_array(as.matrix(get.adjacency(adjacency,names=TRUE,sparse=FALSE))),file=paste(PATH,'network_',(i-1)*1000+j,sep=''))
									}else
									{
										write.table(get.adjacency(adjacency,sparse=FALSE,names=TRUE),file=paste(PATH,'network_',(i-1)*1000+j,sep=''),append=F)
									}
								}else
								{
									if(write.sparse)
									{
										write_stm_CLUTO(as.simple_sparse_array(as.matrix(adjacency)),file=paste(PATH,'network_',(i-1)*1000+j,sep=''))
									}else
									{
										write.table(adjacency,file=paste(PATH,'network_',(i-1)*1000+j,sep=''),append=FALSE)
									}
								}

    					}	

			}




	

}




birewire.visual.monitoring.bipartite<-function(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,exact=FALSE,n.networks=100,perplexity=15,sequence=c(1,5,100,"n"),ncol=2,nrow=length(sequence)/ncol,display=TRUE)
{


if(display)
{
	par(mfrow=c(nrow,ncol))
	par(pty="s")
}
	
dist=list()
tsne=list()
ii=1
for( i in sequence)
{
	print(paste("K=",i))
	data_tmp=data
	tot=list(data)
	m=matrix(nrow=n.networks,ncol=n.networks,0)
	for(j in 2:n.networks)
		{
			#print("entro")
			data_tmp=birewire.rewire.bipartite(incidence=data_tmp,  max.iter=i, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
			#print("esco")
			#print(data_tmp)
			tot[[j]]=data_tmp
			for(k in 1:(j-1))
				{
					#print("grafo")
					#print(tot[[k]])
					m[k,j]=m[j,k]=1-birewire.similarity(tot[[k]],tot[[j]])
				}
				

		}
	dist[[ii]]=m	
	tmp=try(tsne(m,whiten=F,perplexity=perplexity))
	if(!is.double(tmp))
		return(list(dist=list(),tsne=list()))
	tsne[[ii]]=tmp
	#tsne[[ii]]=cmdscale(m,eig=TRUE, k=2)$points
	if(display)
		{
			plot(tsne[[ii]],col=colorRampPalette(c("blue", "red"))( n.networks),pch=16,xlab='A.U.',ylab='A.U.',main=paste('k=',i))
			text(x=tsne[[ii]][1,1],y=tsne[[ii]][1,2],label='start')
		}	
	ii=ii+1


}

return(list(dist=dist,tsne=tsne))



}

birewire.visual.monitoring.undirected<-function(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,exact=FALSE,n.networks=100,perplexity=15,sequence=c(1,5,100,"n"),ncol=2,nrow=length(sequence)/ncol,display=TRUE)
{


if(display)
{
	par(mfrow=c(nrow,ncol))
	par(pty="s")
}
dist=list()
tsne=list()
ii=1
for( i in sequence)
{
	print(paste("K=",i))
	data_tmp=data
	tot=list(data)
	m=matrix(nrow=n.networks,ncol=n.networks,0)
	for(j in 2:n.networks)
		{
			data_tmp=birewire.rewire.undirected(data_tmp,  max.iter=i, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
			tot[[j]]=data_tmp
			for(k in 1:(j-1))
				m[k,j]=m[j,k]=1-birewire.similarity(tot[[k]],tot[[j]])
				

		}
	dist[[ii]]=m	
	tmp=try(tsne(m,whiten=F,perplexity=perplexity))
	if(!is.double(tmp))
		return(list(dist=list(),tsne=list()))
	tsne[[ii]]=tmp
	#tsne[[ii]]=cmdscale(m,eig=TRUE, k=2)$points
	if(display)
		{
			plot(tsne[[ii]],col=colorRampPalette(c("blue", "red"))( n.networks),pch=16,xlab='A.U.',ylab='A.U.',main=paste('k=',i))
			text(x=tsne[[ii]][1,1],y=tsne[[ii]][1,2],label='start')
		}	
	ii=ii+1


}

return(list(dist=dist,tsne=tsne))



}




















##################DSG STUFF#############################
birewire.sampler.dsg<-function(dsg,K,path,delimitators=list(negative='-',positive='+'),exact=FALSE,verbose=TRUE, max.iter.pos='n',max.iter.neg='n', accuracy=0.00005,MAXITER_MUL=10)
	{

  
		if(!is.list(dsg)  )
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

birewire.rewire.dsg<-function(dsg,exact=FALSE,verbose=1,max.iter.pos='n',max.iter.neg='n',accuracy=0.00005,MAXITER_MUL=10,path=NULL,delimitators=list(positive='+',negative= '-'))
{


		if(!is.list(dsg) )
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
if(is.null(table))
	return(NULL)

source=rownames(table)
target=colnames(table)
index=which(table>0,arr.ind=TRUE)
df=data.frame(source=source[index[,1]],sign=sign,target=target[index[,2]])
return(df)

}

##from a sif file, the routine generates the negative and positive incidence matrix
birewire.induced.bipartite<-function(g,delimitators=list(negative='-',positive='+'),sparse=FALSE)
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
	if(!sparse)
	{
		positive=as.data.frame.matrix(table(g_p))
		negative=as.data.frame.matrix(table(g_n))
		dsg[['positive']]=simplify.table(positive)
		dsg[['negative']]=simplify.table(negative)
	}else
	{
			g_p[,2]=paste(g_p[,2],"@@t",sep='')
			g_n[,2]=paste(g_n[,2],"@@t",sep='')
			dsg[['positive']]=graph.edgelist(as.matrix(g_p),directed=TRUE)
			dsg[['negative']]=graph.edgelist(as.matrix(g_n),directed=TRUE)
			V(dsg[['positive']])$type=0
			V(dsg[['positive']])$type[which(unlist(lapply(strsplit(V(dsg[['positive']])$name,"@@"),length))==2)]=1
			V(dsg[['negative']])$type=0
			V(dsg[['negative']])$type[which(unlist(lapply(strsplit(V(dsg[['negative']])$name,"@@"),length))==2)]=1
			

	}
return(dsg)

}	


##inverse of the function above
birewire.build.dsg<-function(dsg,delimitators=list(negative='-',positive='+'))
{


	if(!is.list(dsg))
			    {
			    	stop("The input must be a dsg object (see References) \n")
			    	return(0)
			    }
	positive=dsg[['positive']]
	negative=dsg[['negative']]
	if(!is.igraph(positive))
	{
		g_p=get.data.frame.from.incidence(positive,delimitators[['positive']])
		g_n=get.data.frame.from.incidence(negative,delimitators[['negative']])
		
		}else
		{
			V(positive)$name=unlist(lapply(strsplit(V(positive)$name,"@@"),function(x){return(x[1])}))
			g_n=NULL
			if(!is.null(negative))
				{
					V(negative)$name=unlist(lapply(strsplit(V(negative)$name,"@@"),function(x){return(x[1])}))
					g_n=get.edgelist(negative,names=TRUE)
					g_n=cbind(g_n,delimitators[['negative']])
					g_n=g_n[,c(1,3,2)]
				}
			g_p=get.edgelist(positive,names=TRUE)
			g_p=cbind(g_p,delimitators[['positive']])
			g_p=g_p[,c(1,3,2)]

		
		}
	g=rbind(g_p,g_n)
	return(g)
}

birewire.load.dsg<-function(path)
	{
		

		return(unique(read.table(path,stringsAsFactors=F)))


	}
birewire.save.dsg<-function(g,file)
	{

		
		write.table(g,file,col.names=FALSE,row.names=FALSE,quote=FALSE)

		


	}

##jaccard index for dsg
	birewire.similarity.dsg<-function(m1,m2)
{


if(is.igraph(m1[["positive"]]))
	{
		e.p=length(E(m1[["positive"]]))
		x.p=length(E(graph.intersection(as.undirected(m1[["positive"]]),as.undirected(m2[["positive"]]))))
		e.n=0
		x.n=0
		if(!is.null(m1[["negative"]]))
			{
				e.n=length(E(m1[["negative"]]))
			 	x.n=length(E(graph.intersection(as.undirected(m1[["negative"]]),as.undirected(m2[["negative"]]))))
			 }
		return((x.p+x.n)/(2*e.p+2*e.n-x.p-x.n))
	}else
	{
		x=sum(m1[['positive']]*m2[['positive']]) +sum(m1[['negative']]*m2[['negative']] )
		e=sum(m1[['positive']])+sum(m1[['negative']])
  		return( x/(2*e-x))

	}
  
}


birewire.analysis.dsg<-function(dsg, step=10, max.iter.pos='n',max.iter.neg='n',accuracy=0.00005,verbose=TRUE,MAXITER_MUL=10,exact=FALSE,n.networks=50,display=TRUE)
{
  
		if(!is.list(dsg) )
			    {
			    	stop("The input must be a dsg object (see References) \n")
			    	return(0)
			    }

	incidence_pos=dsg[["positive"]]
	incidence_neg=dsg[["negative"]]
	incidence_pos=birewire.analysis.bipartite(incidence=incidence_pos,  max.iter=max.iter.pos, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact,n.networks=n.networks,display=FALSE,step=step)
    incidence_neg=birewire.analysis.bipartite(incidence=incidence_neg,  max.iter=max.iter.neg, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact,n.networks=n.networks,display=FALSE,step=step)
	

    if(!is.null(incidence_pos$data) & !is.null(incidence_neg$data))
    {
	if(ncol(incidence_pos$data)<ncol(incidence_neg$data))
			{
				mag=incidence_neg
				min=incidence_pos
				mag_is_pos="negative"
				min_is_pos="positive"

				}else
					{
						mag=incidence_pos
						min=incidence_neg
						mag_is_pos="positive"
						min_is_pos="negative"
					}
					}else
					{
						if(is.null(incidence_pos$data))
							 {
								mag=incidence_neg
								min=incidence_pos
								mag_is_pos="negative"
								min_is_pos="positive"

							 }else
							 {
							 	mag=incidence_pos
								min=incidence_neg
								mag_is_pos="positive"
								min_is_pos="negative"
							 }


					}
	if(display)
	{
		try(dev.off())
		#mean=colMeans(mag$data)
		#std=apply(mag$data,2,sd)
		#sup=mean+ qt(.975,nrow(mag$data)-1)*std/sqrt(nrow(mag$data))
		#inf=mean- qt(.975,nrow(mag$data)-1)*std/sqrt(nrow(mag$data))
		#par(mfrow=c(2,1))
		#x=seq(1,length.out=length(mean))
		#plot(step*x,mean,type= 'n',col='blue',lwd=2,main="Jaccard index (JI) over time",xlab="Switching steps",ylab='Jaccard Index',ylim=c(min(mag$data,min$data),max(mag$data,min$data) ))
		#polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey80', border = NA)
		#lines(step*x,mean,col='blue',lwd=2)
		#mean=colMeans(min$data)
		#std=apply(min$data,2,sd)
		#sup=mean+ qt(.975,nrow(min$data)-1)*std/sqrt(nrow(min$data))
		##inf=mean- qt(.975,nrow(min$data)-1)*std/sqrt(nrow(min$data))
		#x=seq(1,length.out=length(mean))
		#polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey60', border = NA)
		#lines(step*x,mean,col='green',lwd=2)
		#abline(v=mag$N,col= 'red')
		#abline(v=min$N,col= 'black')
		#legend("topright",ncol=2,cex=0.8,
		#		col=c('blue','grey80','green','grey60','red','black'),
		#		lwd=c(2,10,2,10,1,1),
		#		legend=c( paste("Mean JI",mag_is_pos), paste("C.I.",mag_is_pos),
		#				  paste("Mean JI",min_is_pos), paste("C.I.",min_is_pos)
		#			,paste("Bound",mag_is_pos), paste("Bound",min_is_pos)))




		mean=colMeans(mag$data)
		std=apply(mag$data,2,sd)
		sup=mean+ qt(.975,nrow(mag$data)-1)*std/sqrt(nrow(mag$data))
		inf=mean- qt(.975,nrow(mag$data)-1)*std/sqrt(nrow(mag$data))
		x=seq(1,length.out=length(mean))
		plot(step*x,mean,type= 'n',col='blue',lwd=2,main="Jaccard index (JI) over time (log-log scale)",log='xy',xlab="Switching steps",ylab='Jaccard Index',ylim=c(min(mag$data[mag$data>0],min$data[min$data>0]),max(mag$data,min$data) ))
		polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey80', border = NA)
		lines(step*x,mean,col='blue',lwd=2)
		if(!is.null(min$data))
		{
		mean=colMeans(min$data)
		std=apply(min$data,2,sd)
		sup=mean+ qt(.975,nrow(min$data)-1)*std/sqrt(nrow(min$data))
		inf=mean- qt(.975,nrow(min$data)-1)*std/sqrt(nrow(min$data))
		x=seq(1,length.out=length(mean))
		polygon(c(rev(step*x),step*x),c(rev(sup),inf), col = 'grey60', border = NA)
		lines(step*x,mean,col='green',lwd=2)
		abline(v=mag$N,col= 'red')
		abline(v=min$N,col= 'black')
		legend("bottomleft",ncol=2,cex=0.8,
				col=c('blue','grey80','green','grey60','red','black'),
				lwd=c(2,10,2,10,1,1),
				legend=c( paste("Mean JI",mag_is_pos), paste("C.I.",mag_is_pos),
						  paste("Mean JI",min_is_pos), paste("C.I.",min_is_pos)
					,paste("Bound",mag_is_pos), paste("Bound",min_is_pos)))
		}else
		{
			abline(v=mag$N,col= 'red')
		legend("bottomleft",ncol=2,cex=0.8,
				col=c('blue','grey80','black'),
				lwd=c(2,10,2,10,1,1),
				legend=c( paste("Mean JI",mag_is_pos), paste("C.I.",mag_is_pos),
					paste("Bound",mag_is_pos)))
		}
		

	}
	return( list(N=list(positive=incidence_pos$N,negative=incidence_neg$N),data=list(positive=incidence_pos$data,negative=incidence_neg$data)))
}


birewire.visual.monitoring.dsg<-function(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,exact=FALSE,n.networks=100,perplexity=15,
	sequence.pos=c(1,5,100,"n"),
	sequence.neg=c(1,5,100,"n"),ncol=2,nrow=length(sequence.pos)/ncol,display=TRUE)
{
if(length(sequence.pos)!=length(sequence.neg))
{
	stop("The two sequence to test must have the same length \n")

}

if(display)
{
	par(mfrow=c(nrow,ncol))
	par(pty="s")
}
dist=list()
tsne=list()
ii=1

for( i in 1:length(sequence.pos))
{
	print(paste("K.pos=",sequence.pos[i],", K.neg=", sequence.neg[i]))
	data_tmp=data
	tot=list(data)
	m=matrix(nrow=n.networks,ncol=n.networks,0)
	for(j in 2:n.networks)
		{
			data_tmp=birewire.rewire.dsg(data_tmp,  max.iter.pos=sequence.pos[i],max.iter.neg=sequence.neg[i], accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
			tot[[j]]=data_tmp
			for(k in 1:(j-1))
				m[k,j]=m[j,k]=1-birewire.similarity.dsg(tot[[k]],tot[[j]])
				

		}
	dist[[ii]]=m	
	tmp=try(tsne(m,whiten=F,perplexity=perplexity))
	if(!is.double(tmp))
		return(list(dist=list(),tsne=list()))
	tsne[[ii]]=tmp
	#tsne[[ii]]=cmdscale(m,eig=TRUE, k=2)$points
	if(display)
		{
			plot(tsne[[ii]],col=colorRampPalette(c("blue", "red"))( n.networks),pch=16,xlab='A.U.',ylab='A.U.',main=paste("K.pos=",sequence.pos[i],", K.neg=", sequence.neg[i]))
			text(x=tsne[[ii]][1,1],y=tsne[[ii]][1,2],label='start')
		}	
	ii=ii+1


}

return(list(dist=dist,tsne=tsne))



}

birewire.slum.to.sparseMatrix<-function(simple_triplet_matrix_sparse) {
  retval <-  sparseMatrix(i=as.numeric(simple_triplet_matrix_sparse$i),
                          j=as.numeric(simple_triplet_matrix_sparse$j),
                          x=as.numeric(as.character(simple_triplet_matrix_sparse$v)),
                          dims=c(simple_triplet_matrix_sparse$nrow, 
                                 simple_triplet_matrix_sparse$ncol),
                          dimnames = dimnames(simple_triplet_matrix_sparse),
                          giveCsparse = TRUE)
  return(retval)
}