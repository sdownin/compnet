##
##  CHECKING k-cycles measure
##

library(igraph)
library(intergraph)
library(ergm)
library(texreg)

g1 <- graph(c(1,2, 1,3, 2,3, 3,4, 4,5),
            directed=F)
  
g2 <- graph(c(1,2, 2,3, 3,4, 4,5),
            directed=F)


n1 <- asNetwork(g1)
n2 <- asNetwork(g2)

mod1 <- ergm(n1 ~ edges + cycle(3) + cycle(4))
mod2 <- ergm(n2 ~ edges + cycle(3) + cycle(4))


mcmc.diagnostics(mod1)
mcmc.diagnostics(mod2)

screenreg(list(mod2))


sna::kcycle.census(n1, maxlen = 4)
sna::kcycle.census(n2, maxlen = 4)

ergm::summary.statistics.ergm(mod2)

summary(n1 ~ edges + cycle(3) + cycle(4), directed=F)

ergmMPLE(n1 ~ edges + cycle(3) + cycle(4), output="array")
ergmMPLE(n2 ~ edges + cycle(3) + cycle(4), output="array")
