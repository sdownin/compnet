##
##
##

library(igraph)
library(intergraph)
library(ergm)
library(texreg)
library(fastcluster)

data_dir <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\firm_nets_rnr"
nets <- readRDS(file.path(data_dir, "qualtrics_d3.rds"))

net <- nets[[length(nets)]]
g <- asIgraph(net)

# M <- net[,]
# cl <- hclust(M)

hc <- fastcluster::hclust(as.dist(net[,]), "average")
plot(hc)

hm <- heatmap(net[,])
plot(hm)


##--------------------
##
##--------------------
## 3 cliques of 4 firms
cliques <- c(1,2, 1,3, 1,4, 2,3, 2,4, 3,4,
             5,6, 5,7, 5,8, 6,7, 6,8, 7,8, 
             9,10, 9,11, 9,12, 10,11, 10,12, 11,12)
## NO complexity
gno <- graph(cliques,
            directed=F)

## LOW complexity
glo <- graph(c(cliques, 1,6, 4,8,  7,11, 5,10,  2,12),
            directed=F)

## HIGH complexity
ghi <- graph(c(cliques, 
               1,6, 1,7, 2,5, 4,5, 4,8,  
               5,9, 5,10, 6,10, 6,9, 7,11, 8,11, 8,12,
               1,9, 2,12, 3,10, 3,11),
            directed=F)

ano <- as_adjacency_matrix(gno, sparse = F)
alo <- as_adjacency_matrix(glo, sparse = F)
ahi <- as_adjacency_matrix(ghi, sparse = F)

hmno <- heatmap(ano)
hmlo <- heatmap(alo)
hmhi <- heatmap(ahi)
  
cols <- c(rep(rgb(0.8,0.1,0.1,0.5),4),
          rep(rgb(0.1,0.8,0.1,0.5),4),
          rep(rgb(0.1,0.1,0.8,0.5),4))
par(mfrow=c(2,2), mar=c(.1,.1,.1,.1))
plot(gno, vertex.color=cols)
plot(glo, vertex.color=cols)
plot(ghi, vertex.color=cols)
  
nno <- asNetwork(gno)
nlo <- asNetwork(glo)
nhi <- asNetwork(ghi)
  
dno <- ergmMPLE(nno ~ edges + cycle(3) + cycle(4), output="array")
dlo <- ergmMPLE(nlo ~ edges + cycle(3) + cycle(4), output="array")
dhi <- ergmMPLE(nhi ~ edges + cycle(3) + cycle(4), output="array")

df <- data.frame(
  cycle3=sapply(list(dno,dlo,dhi),function(x)summary(x[,,2]))
  cycle4=c()
)

df3c <- sapply(list(dno,dlo,dhi),function(x){
    mat <- x$predictor[,,2] ## 4cycle is 2nd predictor
    return(summary(mat[upper.tri(mat,diag = F)]))
  })

df4c <- sapply(list(dno,dlo,dhi),function(x){
  mat <- x$predictor[,,3]  ## 4cycle is 3rd predictor
  return(summary(mat[upper.tri(mat,diag = F)]))
})

df3c
df4c


fno <- ergmMPLE(nno ~ edges + cycle(3) + cycle(4), output="fit")
flo <- ergmMPLE(nlo ~ edges + cycle(3) + cycle(4), output="fit")
fhi <- ergmMPLE(nhi ~ edges + cycle(3) + cycle(4), output="fit")
screenreg(list(fno,flo,fhi))



par(mfrow=c(3,3),mar=c(.1,.1,.1,.1))

## TEST SIMULATION FUNCTION 
testSim <- function(coef, formula=nno ~ edges + cycle(3) + cycle(4), nsim=9) 
{
  ##SIMULATE
  sim1 <- simulate(formula, nsim=nsim, coef = coef)
  ## PLOT
  plot.title <- sprintf("test_cycle_3c_%s_4c_%s.png", 
                        str_replace_all(coef[2], "\\.","_"), 
                        str_replace_all(coef[3], "\\.","_"))
  png(plot.title, width = 8, height = 8, units = 'in', res = 200)
    par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
    sapply(sim1, function(net){
      plot(asIgraph(net), layout=layout.fruchterman.reingold)
    })
    mtext(sprintf("3c=%.2f; 4c=%.2f", coef[2], coef[3]), outer = TRUE, cex = 1.5)
  dev.off()
  ## SUMMARY
  summary(sim1)
  return(sim1)
}

setwd('C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2')

## HIGH 3-cycles | HIGH 4 cycles
hh <- testSim(c(-1.8, 1, .5)) 
## HIGH 3-cycles | LOW 4 cycles
hl <- testSim(c(-1.8, 1, .1)) 
## HIGH 3-cycles | NEGATIVE 4 cycles
hn <- testSim(c(-1.8, 1, -.1)) 
## LOW 3-cycles | HIGH 4 cycles
lh <- testSim(c(-1.8, .2, .5)) 
## LOW 3-cycles | LOW 4 cycles
ll <- testSim(c(-1.8, .2, .1)) 
## LOW 3-cycles | NEGATIVE 4 cycles
ln <- testSim(c(-1.8, .2, -.1)) 


zz <- testSim(c(-1.8, -10, 1), 
              formula=nno ~ edges + cycle(3) + cycle(4)) 




xx <- testSim(c(-1.8, -1, 2)) 









