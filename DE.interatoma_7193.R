
##Step 1: load RedeR package
library(RedeR)
library(igraph)
#data(hs.inter)
#load(file="Data/SH.limmaREF.RData")

dt_12w_12w <- read.csv(file = "/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_12wObeseX12wLean_enriched.csv",header = T)
rownames(dt_12w_12w) <- dt_12w_12w$probes
dt_12w_12w[,1] <- NULL
dt_12w_12w <- dt_12w_12w[!dt_12w_12w$ensembl_peptide_id=="",]
colnames(dt_12w_12w)[4] <- "symbol"

dt_network <- read.table("/home/eulle/Documents/Mestrado/MicroArray/Data/10116.protein.links.v11.0.txt",header = T , stringsAsFactors = F)
dt_network_ <- dt_network[dt_network$combined_score>700,]


dt_network_["id"]<-c(sapply(strsplit(as.character(dt_network_$protein1), split = "\\.") , "[" ,2))
dt_network_["id2"]<-c(sapply(strsplit(as.character(dt_network_$protein2), split = "\\.") , "[" ,2))
dt_network_<-dt_network_[,-c(1,2)]
colnames(dt_network_)<-c("combined_score","protein1","protein2")
rownames(dt_network_)<-NULL
dt_network_<-dt_network_[,c(2,3,1)]

df.g <- graph_from_data_frame(dt_network_,directed = F)

##Step 2: get the subgraphs
gt12w_12w <- subg(g=df.g, dat=dt_12w_12w, refcol=7, maincomp=TRUE, connected=TRUE)
gt12w_12w <- att.setv(g=gt12w_12w, from="symbol", to="nodeAlias")
gt12w_12w <- att.setv(g=gt12w_12w, from="pvalue", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)
E(gt12w_12w)$edgeWidth<-3
V(gt12w_12w)$nodeFontSize<-18

##Step 4: build the server port
rdp <- RedPort ()
calld(rdp)

##Step 3: plot the grafhs
n0 <- addGraph(rdp, gt12w_12w, gcoord=c(15,25), gscale=15, isNest=TRUE, theme='tm1', zoom=25)


##Step 4: nest the grafhs
nestNodes(rdp, nodes=V(gt12w_12w)$name,parent=n0, theme='tm2')


mergeOutEdges(rdp)
relax(rdp,20,100,100,100,10,75,5,20,ps=TRUE)
