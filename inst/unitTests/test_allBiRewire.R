library (BiRewire)
run.tests = function ()
{
test_birewire.analysis()
test_birewire.bipartite()
test_birewire.bipartite.incidence()
test_birewire.bipartite.sparse()
test_birewire.undirected ()
test_birewire.dsg()
}

test_birewire.analysis <- function() 
{
g <- simplify(graph.bipartite( rep(0:1,length=100),
                               c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
                                 seq(1,100,17),seq(1,100,19),seq(1,100,23),100
                               )))
 m<-as.matrix(get.incidence(graph=g))
step=1
max=100*length(E(g))
checkTrue(is.numeric(birewire.analysis(m,step,max)$N))
checkTrue(is.numeric(birewire.analysis(m,step,"n")$N))
}

test_birewire.bipartite <- function(){
g <- simplify(graph.bipartite( rep(0:1,length=100),
                               c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
                                 seq(1,100,17),seq(1,100,19),seq(1,100,23),100
                               )))
 m<-as.matrix(get.incidence(graph=g))
max=100*length(E(g))
step=1
maxiter=birewire.analysis(m,step,max)$N
#checkException(birewire.rewire.bipartite(m,maxiter))
 checkTrue(is.data.frame(birewire.rewire.bipartite(m,maxiter))|is.matrix(birewire.rewire.bipartite(m,maxiter)) )

}

test_birewire.bipartite.incidence <- function(){
g <- simplify(graph.bipartite( rep(0:1,length=100),
                               c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
                                 seq(1,100,17),seq(1,100,19),seq(1,100,23),100
                               )))
 m<-as.matrix(get.incidence(graph=g))
 checkTrue(is.igraph( birewire.bipartite.from.incidence(m,T,F)))
}
test_birewire.bipartite.sparse <- function(){
g <- simplify(graph.bipartite( rep(0:1,length=100),
                               c(c(1:100),seq(1,100,3),seq(1,100,7),100,seq(1,100,13),
                                 seq(1,100,17),seq(1,100,19),seq(1,100,23),100
                               )))
#checkException(birewire.rewire.sparse.bipartite(g))
 checkTrue(is.igraph(birewire.rewire.bipartite(g)))

}
test_birewire.undirected <- function(){
g<-erdos.renyi.game(directed=F,loops=F,n=100,p.or.m=0.2)
  m<-get.adjacency(g,sparse=FALSE)
step=10
max=100*length(E(g))
 checkTrue(is.numeric(birewire.analysis.undirected(m,step,max,verbose=TRUE)$N))
#checkException(birewire.rewire(m,verbose=TRUE))
 checkTrue(is.data.frame(birewire.rewire(m,verbose=TRUE))|is.matrix(birewire.rewire(m,verbose=TRUE)) )
 checkTrue(is.igraph(birewire.rewire(g,verbose=TRUE)))



}

test_birewire.dsg<-function()
{
##g=birewire.load.dsg("data/test.sif")
data(test_dsg)
dsg=birewire.induced.bipartite(test_dsg,delimitators=list(negative='-',positive='+'))
tmp= birewire.rewire.dsg(dsg,verbose=FALSE)
dsg2=birewire.build.dsg(tmp,delimitators=list(negative='-',positive='+'))
##birewire.save.dsg(dsg2,"test2.sif")
checkTrue(birewire.similarity.dsg(dsg,tmp)<=1)

}

