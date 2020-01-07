### Mn minerals data prep for network diagrams ###

library(maps)
library(rjson)
library(readr)
library(plyr)


## read data
# Mn.ages.and.localities <- read.csv("C:/Users/siu854098256/Dropbox/Dan's stuff/Data/Mn minerals/Mn ages and localities.csv",stringsAsFactors = F)
Mn_ages_and_localities <- read_csv("~/Downloads/Mn/Mn ages and localities.csv")
Mn_mineral_groups <- as.data.frame(read_csv("~/Downloads/Mn/Mn_mineral_groups.csv"))
MnMinerals3 <- read_csv("~/Downloads/Mn/MnMinerals3.csv")

## clean data
Mn_localities_minerals <- Mn_ages_and_localities

## separate minerals per locality
# Mn.mins <- strsplit(Mn_localities_minerals$Query.Minerals.at.Locality...,',')
Mn.mins <- strsplit(Mn_localities_minerals$`Query Minerals at Locality`,',')

# Mn_localities_minerals <- data.frame(locality=rep(Mn_localities_minerals$Locality.containing.Mineral..., sapply(Mn.mins, FUN=length)), mineral=unlist(Mn.mins),stringsAsFactors = F)
Mn_localities_minerals <- data.frame(locality=rep(Mn_localities_minerals$`Locality containing Mineral`, sapply(Mn.mins, FUN=length)), mineral=unlist(Mn.mins), age=rep(Mn_ages_and_localities$`Max Age (Ma)`, sapply(Mn.mins, FUN=length)),stringsAsFactors = F)

# remove leading & trailing spaces from mineral names
Mn_localities_minerals$mineral <- sapply(Mn_localities_minerals$mineral,function(x) {sub("^\\s+|^\\s+$","",x)})


## get granular localities

# country <- sapply(strsplit(Mn_localities_minerals$locality_longname,","),function(x) return(gsub("^\\s+|\\s+$", "", x[length(x)])))
# province <- sapply(strsplit(Mn_localities_minerals$locality_longname,","),function(x) return(gsub("^\\s+|\\s+$", "", x[ifelse(length(x)>1,length(x)-1,length(x))])))
# sub.province <- sapply(strsplit(Mn_localities_minerals$locality_longname,","),function(x) return(gsub("^\\s+|\\s+$", "", x[ifelse(length(x)>2,length(x)-2,ifelse(length(x)>1,length(x)-1,length(x)))])))
# 
# Mn_localities_minerals <- cbind(Mn_localities_minerals,as.data.frame(sub.province))
# Mn_localities_minerals <- cbind(Mn_localities_minerals,as.data.frame(province))
# Mn_localities_minerals <- cbind(Mn_localities_minerals,as.data.frame(country))
# 
# Mn_localities_minerals$sub.province <- as.character(Mn_localities_minerals$sub.province)
# Mn_localities_minerals$province <- as.character(Mn_localities_minerals$province)
# Mn_localities_minerals$country <- as.character(Mn_localities_minerals$country)
# 
# lvl1 <- ifelse(Mn_localities_minerals$country%in%c("USA","Canada","Australia","China","Russia"),Mn_localities_minerals$province,Mn_localities_minerals$country)
# lvl2 <- ifelse(Mn_localities_minerals$country%in%c("USA","Canada","Australia","China","Russia"),Mn_localities_minerals$sub.province,Mn_localities_minerals$province)


## take a copy of the dataset and create binary locality/mineral matrix
Mn.loc.min <- Mn_localities_minerals
Mn.loc.min$locality <- as.character(Mn.loc.min$locality)
Mn.loc.min.mat <- as.data.frame.matrix(table(Mn.loc.min[,1:2]))

## filter out countries from localities
rownames(Mn.loc.min.mat[which(rownames(Mn.loc.min.mat)%in%world.cities$country.etc),])
# Mn.loc.min.mat <- Mn.loc.min.mat[-which(row.names(Mn.loc.min.mat)%in%world.cities$country.etc),]

## generate network structure

# create a list of column names (mineral names)
cols = colnames(Mn.loc.min.mat)
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
apply(Mn.loc.min.mat,1,computeCoexistence)

# create list of values of 1 with length = source/target nodes
value <- rep(1,length(src))

# create links (edges) object by combining source, target, and value lists
Mn.links <- as.data.frame(cbind(src,tar,value))

# rename columns
colnames(Mn.links)[1] <- "source"
colnames(Mn.links)[2] <- "target"

# set columns data types
Mn.links$source <- as.character(Mn.links$source)
Mn.links$target <- as.character(Mn.links$target)
Mn.links$value <- as.numeric(as.character(Mn.links$value))

# aggregate values on unique source/target combinations such that the "value" column contains the instances of coexistence
Mn.links <- aggregate(.~source+target,data=Mn.links,FUN=sum)

## compute mineral counts across localities

# create table of mineral/locality
Mn.min.loc <- as.data.frame(cbind(Mn_localities_minerals[,2],Mn_localities_minerals[,1]),stringsAsFactors = F)

# create binary matrix of mineral/locality
Mn.min.loc.mat <- as.data.frame.matrix(table(Mn.min.loc))

Mn.counts <- rowSums(Mn.min.loc.mat)

Mn.counts1 <- Mn.counts[Mn.links$source]
Mn.counts2 <- Mn.counts[Mn.links$target]

Mn.links <- cbind(Mn.links,Mn.counts1,Mn.counts2)

metric <- c()
apply(Mn.links,1,function(x) {metric <<- c(metric,ifelse(as.numeric(x[3])!=min(as.numeric(x[4]),as.numeric(x[5])),(1 - (as.numeric(x[3])/min(as.numeric(x[4]),as.numeric(x[5])))),0.1)) })

Mn.links <- cbind(Mn.links,metric)

Mn.links <- Mn.links[,-c(3,4,5)]

colnames(Mn.links)[3] <- "value"

# create nodes table
Mn.nodes <- as.data.frame(cbind(rownames(Mn.min.loc.mat),rowSums(Mn.min.loc.mat)),stringsAsFactors = F)

# set column names
colnames(Mn.nodes)[1] = "id"
colnames(Mn.nodes)[2] = "size"

# set node sizes
Mn.nodes$size <- as.numeric(Mn.nodes$size)
Mn.nodes$size <- log2(Mn.nodes$size)
Mn.nodes$size[Mn.nodes$size==0] <- 1

# set node groups
rownames(Mn_mineral_groups) <- Mn_mineral_groups$mineral
groups <- Mn_mineral_groups$group[Mn.nodes$id]
groups <- Mn_mineral_groups[Mn.nodes$id,"group"]
Mn.nodes <- cbind(Mn.nodes,groups)
colnames(Mn.nodes)[3] <- "group"

# add paragneteic modes
Mn.pmodes <- as.data.frame(MnMinerals3[which(MnMinerals3$`Mineral Name (plain)`%in%Mn.nodes$id),c("Mineral Name (plain)","Paragenetic mode")])
# Mn.nodes$id[which(!Mn.nodes$id%in%MnMinerals3$`Mineral Name (plain)`)]
rownames(Mn.pmodes) <- Mn.pmodes$`Mineral Name (plain)`

pmodes <- Mn.pmodes[Mn.nodes$id,"Paragenetic mode"]

pmodes1 <- sapply(pmodes,function(x){ ifelse(length(strsplit(x,",")[[1]])>1,substr(as.character(max(sapply(strsplit(x,","),function(x){as.numeric(sub('"',"",x))}))),1,1),substr(x,1,1)) })

Mn.nodes <- cbind(Mn.nodes,pmodes1)

colnames(Mn.nodes)[3] <- "Redox"

# add age
Mn.ages <- Mn_localities_minerals[,2:3]

ages <- aggregate(.~mineral,data = Mn.ages,FUN = max)

Mn.nodes <- cbind(Mn.nodes,apply(Mn.nodes,1,function(x){x["age"]<-ifelse(x["id"]%in%ages$mineral,ages$age[which(ages$mineral==x["id"])],NA)}))

colnames(Mn.nodes)[5] <- "age"

Mn.nodes$age <- as.numeric(Mn.nodes$age)

Mn.nodes$Age_Intervals<-cut(x = Mn.nodes$age,breaks = c(0,349,499,799,1199,1599,1999,2999,3999,max(Mn.nodes$age)))

#Mn.Node_Colors<-randomColor(count = nlevels(Mn.nodes$Age_Intervals))

Mn.Node_Colors<- c("red","blue","green","yellow","orange","white","cyan","purple")


Mn.nodes$color <- Mn.Node_Colors[as.numeric(Mn.nodes$Age_Intervals)]
Mn.nodes$color[is.na(Mn.nodes$age)] <- "lime"

cut2()

## less density
Mn.links.07 <- Mn.links[which(Mn.links$value<0.7),]



## write output files
library(jsonlite)
Mn.network.json = paste("{ \"nodes\": ",toJSON(unname(split(Mn.nodes, 1:nrow(Mn.nodes)))),",\"links\": ",toJSON(unname(split(Mn.links, 1:nrow(Mn.links)))),"}")
write(Mn.network.json, "~/PycharmProjects/dtdi-network-diagrams/Mn/network/Mn_network_new.json")
