library(limma)
library(affy)

#Pré-Processamento

# Importar a tabela da plataforma e filtrar pelas colunas de interesse
df_platform <- read.delim("/home/eulle/Documents/Mestrado/MicroArray/Datasets/GPL1355-10794.txt",header = T,stringsAsFactors = F)
df_platform <- df_platform[,c("ID","Gene.Symbol","ENTREZ_GENE_ID")]

# Sondas que possuem mais de uma referência.
symbol <- strsplit(df_platform$Gene.Symbol,split = "///")
entrez <- strsplit(df_platform$ENTREZ_GENE_ID,split = "///")
df_platform_symbol <- data.frame(ID = rep(df_platform$ID, sapply(symbol, length)), Gene.Symbol = unlist(symbol))
colnames(df_platform_symbol) <- c("probes","external_gene_name")

df_platform_entrezid <- data.frame(ID = rep(df_platform$ID, sapply(entrez, length)), ENTREZ_GENE_ID = unlist(entrez))
df_platform[df_platform==""] <- NA

library(biomaRt)

# mart <- useMart("ENSEMBL_MART_ENSEMBL",host="http://uswest.ensembl.org/")
# mart <- useDataset("rnorvegicus_gene_ensembl",mart)
# external<-as.vector(df_platform_symbol[,1])
# annotLookup <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),filters = "external_gene_name",values = external,mart = mart)

ensembl <- useMart("ensembl",dataset = "rnorvegicus_gene_ensembl",host="http://uswest.ensembl.org/")
gene_map <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name","ensembl_peptide_id"),filters = "external_gene_name",values = df_platform_symbol$external_gene_name,mart = ensembl)
# Fiz o merge para pegar a anotação dos emsembl_gene_id e retirei os valores duplicados das sondas. Neste processo alguns genes
# foram descartados pois não tinham nenhuma informação relevante,como o ensemble_gene_id, ensemble_transcript_id...
df_gene_map <- merge(gene_map,df_platform_symbol,by="external_gene_name",all.y = TRUE)
df_gene_map_aduplicated <- df_gene_map[!duplicated(df_gene_map$probes) & !is.na(df_gene_map$ensembl_gene_id) ,]

# Pegar todos os valores repetidos e realizar um cálculo de dispersão com os dados de expressão.
allDup <- function (value) 
{ 
  duplicated(value) | duplicated(value, fromLast = TRUE) 
} 

df_ensembl_duplicated <- df_gene_map_aduplicated[allDup(df_gene_map_aduplicated$ensembl_gene_id),]
rownames(df_ensembl_duplicated) <- NULL
# Realizar a comparação entre os dados duplicados e os dados de expressão.
library(data.table)
library(matrixStats)
celpath = "/home/eulle/Documents/Mestrado/MicroArray/Datasets/GSE7253_RAW/"

#celpath_obese_lean = '/home/eulle/Documents/Mestrado/Microarray/Datasets/GSE7193_RAW/'

dt = ReadAffy(celfile.path = celpath)


#dt_2 = ReadAffy(celfile.path = celpath_obese_lean)

# Converter um AffyBatch para um ExpressionSet object
#Normalização
exp <- rma(dt)

#exp_7193 <- rma(dt_2)

data <- exprs(exp)

data<-as.data.frame(data)
data<-setDT(data,keep.rownames = TRUE)
# probes_ensembl são os valores de sondas que são referenciados pelo mesmo ensembl_gene_id, então é necessário o cálculo de desvio padrão.
data_probes_ensembl <- data[data$rn%in%df_ensembl_duplicated$probes,]
data<-data[!data$rn%in%df_ensembl_duplicated$probes,]

unique_ensemble_duplicated <- unique(df_ensembl_duplicated$ensembl_gene_id)
data_sd <- data.frame()
data_sd_rbind <- data_sd


for(ii in 1:length(unique_ensemble_duplicated)){
  df_subset_ensembl <- df_ensembl_duplicated[df_ensembl_duplicated$ensembl_gene_id==unique_ensemble_duplicated[ii],]
  df_subset_data <- data_probes_ensembl[data_probes_ensembl$rn%in%df_subset_ensembl$probes,]
  
  df_subset_data<- data.frame(df_subset_data[,-1],row.names = df_subset_data$rn)
  #df_subset_data <- df_subset_data[,-1]
  data_sd <- as.data.frame(df_subset_data[which.max(rowSds(as.matrix(df_subset_data))),])
  data_sd_rbind <- rbind(data_sd_rbind,data_sd)
}

temp <- df_gene_map_aduplicated[isUnique(df_gene_map_aduplicated[,2]),]
temp_data <- data[data$rn%in%temp$probes,]

data_sd_rbind<-setDT(data_sd_rbind,keep.rownames = TRUE)
data_rbind_all<-rbind(temp_data,data_sd_rbind)
data_rbind_all<-data.frame(data_rbind_all[,-1],row.names = data_rbind_all$rn)

## Armazenar os dados


#data_7193 <- exprs(exp_7193)

# library(biomaRt)
# #dt_ensembl <- read.table(file = "/home/eulle/Documents/Mestrado/MicroArray/mart_export.txt",header = T,na.strings = "NA" ,sep = ",")
# #mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl" , host = "www.ensembl.org")
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL",host="http://uswest.ensembl.org/")
# mart <- useDataset("rnorvegicus_gene_ensembl",mart)
# #g <- getGene( id = "1370583_s_at",type = "affy_rat230_2", mart = mart)
# annotLookup <- getBM(mart=mart, attributes=c("affy_rat230_2", "ensembl_peptide_id","ensembl_transcript_id","ensembl_gene_id" ,"gene_biotype"), filter="affy_rat230_2", values=rownames(data), uniqueRows=TRUE)

#colnames(data) <- c(rep("control_10w", 3), rep("diabetic_10w", 3), rep("control_20w", 3), rep("diabetic_20w",3))
tps <- c(rep("control_10w", 3), rep("diabetic_10w", 3), rep("control_20w", 3), rep("diabetic_20w",3))

t <- factor(tps)

#ph = dt@phenoData
#list = list.files(celpath,full.names=TRUE)

#Dizer a qual grupo determinada amostra pertence
#treats <- as.factor(colnames(data))

design <- model.matrix(~0+t)

fit <- lmFit(data_rbind_all, design)
# fit$genes$entrez <- entrez
# fit$genes$symbol <- symbol


contrasts <- makeContrasts(tdiabetic_20w-tcontrol_20w, levels = design)

ct.fit <- eBayes(contrasts.fit(fit, contrasts))

res.fit <- decideTests(ct.fit, method = "global", adjust.method = "BH", p.value = 0.05,lfc = 0.05)

degs <- data.frame(LogFC = ct.fit$coefficients, p.value = ct.fit$p.value, degenes = unclass(res.fit),stringsAsFactors = F)

ups <- subset(degs,degs[,3]==1)
downs <- subset(degs,degs[,3]==-1)
degs_up_down <- rbind(downs,ups)
probes_up_down <- rownames(degs_up_down)
genes_up_down <- df_gene_map_aduplicated[df_gene_map_aduplicated$probes%in%probes_up_down,]
degs_up_down <- setDT(degs_up_down,keep.rownames = TRUE)
genes_up_down <- merge(degs_up_down,genes_up_down,by.x ="rn",by.y = "probes",all.y = TRUE)
colnames(genes_up_down)[c(1,2,3,4)] <- c("probes","logFC","pvalue","degs")
#genes_up_down<- genes_up_down[,-c(1,2,3)]
#colnames(genes_up_down["tdiabetic_10w...tcontrol_10w.2"]) <- genes_up_down["degs"]

degs_universe <- degs
probes_universe <- rownames(degs_universe)
genes_universe <- df_gene_map_aduplicated[df_gene_map_aduplicated$probes%in%probes_universe,]
degs_universe <- setDT(degs_universe,keep.rownames = TRUE)
genes_universe <- merge(degs_universe,genes_universe,by.x ="rn",by.y = "probes",all.y = TRUE)
colnames(genes_universe)[c(1,2,3,4)] <- c("probes","logFC","pvalue","degs")
# genes_universe <- genes_univer se[!duplicated(genes_universe$affy_rat230_2),]
# rownames(genes_universe) <- genes_universe$affy_rat230_2

# genes_universe$affy_rat230_2 <- NULL
# geness<- merge(degs_universe,genes_universe,by="row.names",all.x=TRUE)
# #genestest<- na.omit(geness)
# rownames(geness) <- geness$Row.names
# geness$Row.names <- NULL
# #geness<-geness[,-5]
# geness<-geness[complete.cases(geness),]
#geness<-geness[!geness$ensembl_peptide_id=="",]

#degs_up_down <- degs_up_down[complete.cases(degs_up_down),]

# probes_up_down <- rownames(degs_up_down)
# genes_up_down <- df_gene_map_aduplicated[df_gene_map_aduplicated$probes%in%probes_up_down,]
# genes_up_down<-genes_up_down[!duplicated(genes_up_down$affy_rat230_2),]
# rownames(genes_up_down) <- genes_up_down$affy_rat230_2
# genes_up_down$affy_rat230_2 <- NULL
# genes <- merge(degs_up_down,genes_up_down,by="row.names",all.x=TRUE)
# rownames(genes) <- genes$Row.names
# genes$Row.names <- NULL
# genes<-genes[,-5]
# genes<-genes[complete.cases(genes),]
# genes<-genes[!genes$ensembl_peptide_id=="",]

# genes_vector <- as.vector(genes[,5])
# geness_vector <- as.vector(geness[,5])

#colnames(genes) <- c("LogF10w...10w","pValue","degenes","ensembl_peptide_id","symbol")
write.csv(genes_up_down,file="/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_10wdiabeticX10wcontrol.csv",row.names = FALSE)

# Enriquecimento funcional
library(clusterProfiler)
library(rat2302.db)
ego2 <- enrichGO(gene         = genes_up_down$ensembl_gene_id,
                 universe = genes_universe$ensembl_gene_id,
                 OrgDb         = rat2302.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)
enrich_data <- ego2@result
enrich_cutOff <- enrich_data[enrich_data$p.adjust <= 0.05,]

# Pega a coluna do geneID , trata as barras , da um unique e xablau !

enriched_genes <- enrich_cutOff[,8]
enriched_genes <- strsplit(enriched_genes,split = "/")
enriched_genes <- unlist(enriched_genes)
enriched_genes <- unique(enriched_genes)

enriched_genes <- genes_up_down[genes_up_down$ensembl_gene_id%in%enriched_genes]
write.csv(enriched_genes,file="/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_20wdiabeticX20wcontrol_enriched.csv",row.names = FALSE)

# Gerar os gráficos matemáticos 

## Gráfico de pizza para verificar a diferença entre genes enriquecidos e os não enriquecidos

## 10w_10w
slices <- c(length(enriched_genes$external_gene_name),length(genes_up_down$external_gene_name)) 
pct <- round(slices)
lbls <- pct # add percents to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart 10wTreated_10wControl")
legend("topright", c("Enriched","Not enriched"), cex = 0.8, fill = rainbow(length(slices)))

## 20w_20w
slices <- c(length(enriched_genes$external_gene_name),length(genes_up_down$external_gene_name)) 
pct <- round(slices)
lbls <- pct # add percents to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart 20wTreated_20wControl")
legend("topright", c("Enriched","Not enriched"), cex = 0.7, fill = rainbow(length(slices)))

## 12w_12w
slices <- c(length(enriched_genes7193$external_gene_name),length(genes7193_up_down$external_gene_name)) 
pct <- round(slices)
lbls <- pct # add percents to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart 12wTreated_12wControl")
legend("topright", c("Enriched","Not enriched"), cex = 0.7, fill = rainbow(length(slices)))

## Verificar se há mudança na expressão gênica dos genes que são aparecem tanto no estudo 10w_10w quanto no 20w_20w

dt_10w10w = read.csv(file = "/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_10wdiabeticX10wcontrol.csv",stringsAsFactors = F)
dt_20w20w = read.csv(file = "/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_20wdiabeticX20wcontrol.csv",stringsAsFactors = F)

degs_10w_20w <- dt_20w20w[dt_20w20w$external_gene_name%in%dt_10w10w$external_gene_name,]



# Não houve enriqueciemento na condição 10w_tratado VERSUS 10w_controle devido ao baixo número de degs.






# zuera <- enrichKEGG(gene=genes_up_down$external_gene_name, organism = "rno", keyType = "kegg",
#            pvalueCutoff = 0.05, pAdjustMethod = "BH", universe =genes_universe$external_gene_name )

# pca<-prcomp(degs)
# autoplot(pca)

# Refazendo tudo, caralho !

# Carregamento das bibliotecas
# library(limma)
# library(affy)




# data<-data.frame(data[,-1],row.names = data$rn)
# data<-rbind(data,data_sd_rbind)
# data<-setDT(data,keep.rownames = TRUE)




#Retirar as linha enriquecidas com padjust <= 0.05





  #IDs.up <- ups
#IDs.downs <- rownames(downs)

#library("biomaRt")

#aaa <- topTable(ct.fit, number = 200, adjust.method = "BH")

# library(biomaRt)
# dt_ensembl <- read.table(file = "/home/eulle/Documents/Mestrado/MicroArray/mart_export.txt",header = T,na.strings = "NA" ,sep = ",")
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL")
# mart <- useDataset("rnorvegicus_gene_ensembl",mart)
# g <- getGene( id = "1370583_s_at",type = "affy_rat230_2", mart = mart)
# annotLookup <- getBM(mart=mart, attributes=c("affy_rat230_2", "ensembl_gene_id", "gene_biotype", "external_gene_name"), filter="affy_rat230_2", values=rownames(data), uniqueRows=TRUE)

# ################################
# 
# ## Probes up regulated
# 
# probes_ups <- rownames(ups) 
# genes_ups <- annotLookup$affy_rat230_2%in%probes_ups
# genes_ups <- annotLookup[genes_ups,]
# rownames(genes_ups) <- genes_ups$affy_rat230_2
# genes_ups$affy_rat230_2 <- NULL
# 
# 
# rownames(genes_ups) <- genes_ups$Row.names
# genes_ups[,c(1 , 2, 3 , 5 , 6 )] <- NULL
# colnames(genes_ups)[names(genes_ups) == "external_gene_name"] <- "symbol"
# 
# ## Probes down regulated
# 
# probes_downs <- rownames(downs) 
# genes_downs <- annotLookup$affy_rat230_2%in%probes_downs
# genes_downs <- annotLookup[genes_downs,]
# 
# rownames(genes_downs) <- genes_downs$affy_rat230_2
# genes_downs$affy_rat230_2 <- NULL
# downs["1389809_at","symbol"]<- "Pmepa1"
# genes_downs <- downs
# rownames(genes_ups) <- genes_ups$Row.names
# genes_downs[,1] <- NULL
# 
# genes_downs <- genes_downs[complete.cases(genes_downs),]
# degs_up_down <- rbind(genes_downs,genes_ups)
# 
# colnames(degs_up_down) <- c("symbol","LogFC","pValue","degenes")
# 
# write.csv(degs_up_down,file="/home/eulle/Documents/Mestrado/MicroArray/GSE7253-Results/degenes_up_down.csv")
# 
