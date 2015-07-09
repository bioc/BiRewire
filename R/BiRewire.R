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
	for(i in 1:n.networks)
	{
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
	RES=rbind(RES,result$similarity_scores)

	}
	if(display)
	{
		mean=colMeans(RES)
		std=apply(RES,2,sd)
		sup=mean+1.96*std/sqrt(nrow(RES))
		inf=mean-1.96*std/sqrt(nrow(RES))
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
	if(is.igraph(m1))
	{
		if(is.bipartite(m1))
		{
			m1=get.incidence(m1,sparse=F)
			m2=get.incidence(m2,sparse=F)
		}else
		{
			m1=get.adjacency(m1,sparse=F)
			m2=get.adjacency(m2,sparse=F)
		}

	}
	
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
  names=V(g)$label
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
  
  gg<-graph.bipartite(edges=result,types=V(g)$type,directed=FALSE)
  if(!is.null(names))
 	 V(gg)$label=names
  
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
	for( i in 1:n.networks)
	{
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
		RES=rbind(RES,result$similarity_scores)

	}
	if(display)
	{
		mean=colMeans(RES)
		std=apply(RES,2,sd)
		sup=mean+1.96*std/sqrt(nrow(RES))
		inf=mean-1.96*std/sqrt(nrow(RES))
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
						MAXITER_MUL=MAXITER_MUL*max.iter
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
        names=V(graph)$label
        

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
			 		max.iter=ceiling((e/(2*d^3-6*d^2+2*d+2))*log(x=(1-d)/accuracy))
				
 				result<-.Call("R_rewire_sparse", edges,n , n, as.numeric( max.iter),e,verbose,0)

			}

		
        gg<-graph(edges=result+1,directed=FALSE,n=n)
        if(!is.null(names))
        	V(gg)$label=names
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
										write_stm_CLUTO(as.simple_sparse_array(as.matrix(get.adjacency(incidence,sparse=F))),file=paste(PATH,'network_',(i-1)*1000+j,sep=''))
									}else
									{
										write.table(get.adjacency(incidence,sparse=F),file=paste(PATH,'network_',(i-1)*1000+j,sep=''),append=F)
									}
								}else
								{
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




	

}



birewire.visual.monitoring.bipartite<-function(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,exact=FALSE,n.networks=100,perplexity=15,sequence=c(1,5,100,"n"),ncol=2,nrow=length(sequence)/ncol,display=TRUE)
{


if(display)
{
	par(mfrow=c(nrow,ncol))
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
			data_tmp=birewire.rewire.bipartite(data_tmp,  max.iter=i, accuracy=accuracy,verbose=verbose,MAXITER_MUL=MAXITER_MUL,exact=exact)
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

birewire.visual.monitoring.undirected<-function(data,accuracy=0.00005,verbose=FALSE,MAXITER_MUL=10,exact=FALSE,n.networks=100,perplexity=15,sequence=c(1,5,100,"n"),ncol=2,nrow=length(sequence)/ncol,display=TRUE)
{


if(display)
{
	par(mfrow=c(nrow,ncol))
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



