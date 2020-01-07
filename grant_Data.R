#Author : Anirudh Prabhu
#email : prabha2@rpi.edu

grants_ec1 <- read.csv("~/Dropbox/Work/grants_ec1.csv")

summary(grants_ec1)

grants_sub <- data.frame(title = grants_ec1$Title,org = grants_ec1$Organization)
summary(grants_sub)
# if there are duplicates, aggregate them for increased edge weights. 
duplicated(grants_sub)

# create a list of column names (mineral names)
cols = colnames(grants.org.mat)
# create empty lists from source and target nodes
src = c()
tar = c()

# function to iterate over columns finding those minerals that coexist per row (locality)
computeCoexistence <- function(x) {
  for (c1 in 1:length(x)) {
    if (x[c1] > 0) {
      for (c2 in c1:length(x)) {
        if (c1 != c2 && x[c2] > 0) {
          # print(paste(cols[c1],cols[c2]))
          src <<- c(src,cols[c1])
          tar <<- c(tar,cols[c2])
        }
      }
    }
  }
}

### apply custom function to dataset
apply(grants.org.mat,1,computeCoexistence)

# create list of values of 1 with length = source/target nodes
value <- rep(1,length(src))

# create links (edges) object by combining source, target, and value lists
grants.links <- as.data.frame(cbind(src,tar,value))

library(igraph)
g<-graph_from_data_frame(grants.links,directed = F)
plot(g,vertex.label = NA)
grants.org.mat <- as.data.frame.matrix(table(grants_sub[,c(2,1)]))


library(network)
library(intergraph)

g1<-asNetwork(g)
plot(network(g1))

