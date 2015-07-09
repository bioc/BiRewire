library (BiRewire)
run.tests = function ()
{
test_birewire.analysis.bipartite()
test_birewire.bipartite()
test_birewire.bipartite.incidence()
test_birewire.bipartite.sparse()
test_birewire.undirected ()
test_birewire.visual.monitoring()
}

test_birewire.analysis.bipartite <- function() 
{
g <- bipartite.random.game(n1=10,n2=4,p=0.5)
 m<-as.matrix(get.incidence(graph=g))
step=1
max=100*length(E(g))
checkTrue(is.numeric(birewire.analysis.bipartite(m,step,max,display=F,n.networks=3)$N))
checkTrue(is.numeric(birewire.analysis.bipartite(m,step,"n",display=F,n.networks=3)$N))
}

test_birewire.bipartite <- function(){
g <- bipartite.random.game(n1=10,n2=4,p=0.5)

 m<-as.matrix(get.incidence(graph=g))
max=100*length(E(g))
step=1
maxiter=birewire.analysis.bipartite(m,step,max,display=F,n.networks=3)$N
#checkException(birewire.rewire.bipartite(m,maxiter))
 checkTrue(is.data.frame(birewire.rewire.bipartite(m,maxiter))|is.matrix(birewire.rewire.bipartite(m,maxiter)) )

}

test_birewire.bipartite.incidence <- function(){
g <- bipartite.random.game(n1=10,n2=4,p=0.5)
 m<-as.matrix(get.incidence(graph=g))
 checkTrue(is.igraph( birewire.bipartite.from.incidence(m,T)))
}
test_birewire.bipartite.sparse <- function(){
g <- bipartite.random.game(n1=10,n2=4,p=0.5)
#checkException(birewire.rewire.sparse.bipartite(g))
 checkTrue(is.igraph(birewire.rewire.bipartite(g)))

}
test_birewire.undirected <- function(){
g<-erdos.renyi.game(directed=F,loops=F,n=100,p.or.m=0.2)
  m<-get.adjacency(g,sparse=FALSE)
step=10
max=100*length(E(g))
 checkTrue(is.numeric(birewire.analysis.undirected(m,step,max,verbose=TRUE,display=F)$N))
#checkException(birewire.rewire(m,verbose=TRUE))
 checkTrue(is.data.frame(birewire.rewire.undirected (m,verbose=TRUE))|is.matrix(birewire.rewire.undirected(m,verbose=TRUE)) )
 checkTrue(is.igraph(birewire.rewire.undirected (g,verbose=TRUE)))



}

test_birewire.visual.monitoring <- function(){
g <- bipartite.random.game(n1=20,n2=4,p=0.5)
b=birewire.visual.monitoring.bipartite(g,display=F,n.networks=10)
g <- erdos.renyi.game(100,0.1)
b=birewire.visual.monitoring.undirected(g,display=F,n.networks=10)




}
