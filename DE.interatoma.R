##Step 1: load RedeR package
library(RedeR)
library(igraph)

# dt_10w_10w <- read.csv(file = "/home/eulle/Documents/Mestrado/Microarray/MicroArray/results_new_method/genes_20wdiabeticX20wcontrol_enriched.csv",header = T)
# rownames(dt_10w_10w) <- dt_10w_10w$X
# dt_10w_10w[,1] <- NULL

# dt_10w_20w <- read.csv(file = "/home/eulle/Documents/Mestrado/Microarray/MicroArray/Results/genes_10wdiabeticX20wcontrol_p.csv" , header = T)
# rownames(dt_10w_20w) <- dt_10w_20w$X
# dt_10w_20w[,1] <- NULL

# dt_20w_10w <- read.csv(file = "/home/eulle/Documents/Mestrado/Microarray/MicroArray/Results/genes_20wdiabeticX10wcontrol_p.csv" , header = T)
# rownames(dt_20w_10w) <- dt_20w_10w$X
# dt_20w_10w[,1] <- NULL

dt_20w_20w <- read.csv(file = "/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_20wdiabeticX20wcontrol_enriched.csv" , header = T)
rownames(dt_20w_20w) <- dt_20w_20w$probes
dt_20w_20w[,1] <- NULL
dt_20w_20w <- dt_20w_20w[!dt_20w_20w$ensembl_peptide_id=="",]
colnames(dt_20w_20w)[4] <- "symbol"

dt_network <- read.table("/home/eulle/Documents/Mestrado/MicroArray/Data/10116.protein.links.v11.0.txt",header = T , stringsAsFactors = F)
dt_network_ <- dt_network[dt_network$combined_score>700,]


dt_network_["id"]<-c(sapply(strsplit(as.character(dt_network_$protein1), split = "\\.") , "[" ,2))
dt_network_["id2"]<-c(sapply(strsplit(as.character(dt_network_$protein2), split = "\\.") , "[" ,2))
dt_network_<-dt_network_[,-c(1,2)]
colnames(dt_network_)<-c("combined_score","protein1","protein2")
rownames(dt_network_)<-NULL
dt_network_<-dt_network_[,c(2,3,1)]

# library(STRINGdb)
# string_db=STRINGdb$new(version="10",species=10116,score_threshold=400,input_directory="../Data/10116.protein.links.v11.0.txt")
# myDF=string_db$map(dt_10w_10w,"symbol",removeUnmappedRows=TRUE)

# Transformar dtf em igraph
df.g <- graph_from_data_frame(dt_network_,directed = F)

#load(file="Data/SH.limmaREF.RData")

##Step 2: get the subgraphs
# gt10w <- subg(g=df.g, dat=dt_10w_10w, refcol = 4, maincomp=TRUE, connected=TRUE)
# gt10w <- att.setv(g=gt10w, from="symbol", to="nodeAlias")
# gt10w <- att.setv(g=gt10w, from="LogF10w...10w", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)
# E(gt10w)$edgeWidth<-3
# V(gt10w)$nodeFontSize<-18

#   
# gt10w_20w <- subg(g=df.g, dat=dt_10w_20w, refcol=4, maincomp=TRUE, connected=TRUE)
# gt10w_20w <- att.setv(g=gt10w_20w, from="symbol", to="nodeAlias")
# gt10w_20w <- att.setv(g=gt10w_20w, from="LogF10w...20w", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)
# E(gt10w_20w)$edgeWidth<-3
# V(gt10w_20w)$nodeFontSize<-18
# 
# gt20w_10w <- subg(g=df.g, dat=dt_20w_10w, refcol=4, maincomp=TRUE, connected=TRUE)
# gt20w_10w <- att.setv(g=gt20w_10w, from="symbol", to="nodeAlias")
# gt20w_10w <- att.setv(g=gt20w_10w, from="LogF20w...10w", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)
# E(gt20w_10w)$edgeWidth<-3
# V(gt20w_10w)$nodeFontSize<-18

gt20w_20w <- subg(g=df.g, dat=dt_20w_20w, refcol=7, maincomp=TRUE, connected=TRUE)
gt20w_20w <- att.setv(g=gt20w_20w, from="symbol", to="nodeAlias")
gt20w_20w <- att.setv(g=gt20w_20w, from="logFC", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)
E(gt20w_20w)$edgeWidth<-3
V(gt20w_20w)$nodeFontSize<-18

##Step 4: build the server port
rdp <- RedPort ()
calld(rdp)

##Step 3: plot the grafhs
# n0 <- addGraph(rdp, gt10w, gcoord=c(15,25), gscale=15, isNest=TRUE, theme='tm1', zoom=25)
# n1 <- addGraph(rdp, gt10w_20w, gcoord=c(37,20), gscale=30, isNest=TRUE, theme='tm1', zoom=25)
# n2 <- addGraph(rdp, gt20w_10w, gcoord=c(65,20), gscale=30, isNest=TRUE, theme='tm1', zoom=25)
n3 <- addGraph(rdp, gt20w_20w, gcoord=c(20,70), gscale=50, isNest=TRUE, theme='tm1', zoom=25)
# n4 <- addGraph(rdp, gt5d, gcoord=c(60,70), gscale=60, isNest=TRUE, theme='tm1', zoom=25)

##Step 4: nest the grafhs
#nestNodes(rdp, nodes=V(gt10w)$name, parent=n1, theme='tm2')
#nestNodes(rdp, nodes=V(gt10w_20w)$name, parent=n2, theme='tm2')
#nestNodes(rdp, nodes=V(gt20w_10w)$name, parent=n3, theme='tm2')
nestNodes(rdp, nodes=V(gt20w_20w)$name, parent=n3, theme='tm2')

mergeOutEdges(rdp)
relax(rdp,20,100,100,100,10,75,5,20,ps=TRUE)
