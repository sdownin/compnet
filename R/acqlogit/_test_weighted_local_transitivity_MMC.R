library(igraph)

set.seed(1)
gz <- erdos.renyi.game(15, type='gnp', p.or.m =.19, directed = F)
z <- sapply(transitivity(gz, type='barrat'), function(x) 7 + 30*ifelse(is.na(x)|is.nan(x),0,x))
plot(gz, vertex.size=z)

##
#
##
transitivMMC <- function(g) {
  mmc <- transitivity(g, type='barrat')
  return( sapply(mmc, function(x) 5 + 30*ifelse(is.nan(x)|is.na(x), 0, x)) )
}

ratioMMC <- function(g) {
  eids <- incident_edges(g, V(g))
  return(
    sapply(eids, function(l) {
      w <- igraph::get.edge.attribute(g, name = 'weight', index = l)
      n.smc <- length(w[w<=1])
      if (n.smc == 0)
        return(0)
      n.mmc <- length(w[w >1])
      return( n.mmc / n.smc )
    })
  )
}


diffMMC <- function(g) {
  eids <- incident_edges(g, V(g))
  return(
    sapply(eids, function(l) {
      w <- igraph::get.edge.attribute(g, name = 'weight', index = l)
      n.smc <- length(w[w<=1])
      n.mmc <- length(w[w >1])
      return( n.mmc - n.smc )
    })
  )
}


m1 <- matrix(c(0,1,1,1,0,0, 
               1,1,0,1,0,0,
               1,0,0,1,1,1,
               0,0,1,0,0,0,
               0,0,1,0,1,0,
               0,0,0,0,0,1,
               0,0,0,0,0,1), nrow=7, byrow=T)
m2 <- matrix(c(0,1,1,1,1,0, 
               1,1,0,1,0,0,
               1,0,0,1,1,1,
               0,0,1,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,1,
               0,0,0,0,0,1), nrow=7, byrow=T)
gb1 <- igraph::graph_from_incidence_matrix(m1, directed = F)
g1 <- bipartite.projection(gb1, multiplicity = T)$proj1

gb2 <- igraph::graph_from_incidence_matrix(m2, directed = F)
g2 <- bipartite.projection(gb2, multiplicity = T)$proj1

df <- data.frame(v=as.integer(V(g1)), 
                 tran1=transitivMMC(g1),
                 tran1u=transitivity(g1, type='local'),
                 ratio1=ratioMMC(g1),
                 diff1=diffMMC(g1), 
                 tran2=transitivMMC(g2), 
                 tran2u=transitivity(g2, type='local'),
                 ratio2=ratioMMC(g2),
                 diff2=diffMMC(g2),
                 stringsAsFactors = F
                 )
df
cor(df[,-1])

plot(g1, 
     edge.width=E(g1)$weight^3,
     vertex.size=transitivMMC(g1)
     )
plot(g2, 
     edge.width=E(g2)$weight^3,
     vertex.size=transitivMMC(g2)
     )
