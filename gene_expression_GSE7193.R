library(limma)
library(affy)

#data("rat2302cdf")
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

## Armazenar os dados
celpath_7193 = "/home/eulle/Documents/Mestrado/MicroArray/Datasets/GSE7193_RAW/"

#celpath_obese_lean = '/home/eulle/Documents/Mestrado/Microarray/Datasets/GSE7193_RAW/'

dt_7193 = ReadAffy(celfile.path = celpath_7193)

#dt_2 = ReadAffy(celfile.path = celpath_obese_lean)

# Converter um AffyBatch para um ExpressionSet object
#Normalização
exp_7193 <- rma(dt_7193)

#exp_7193 <- rma(dt_2)

data_7193 <- exprs(exp_7193)

data_7193<-as.data.frame(data_7193)
data_7193<-setDT(data_7193,keep.rownames = TRUE)
# probes_ensembl são os valores de sondas que são referenciados pelo mesmo ensembl_gene_id, então é necessário o cálculo de desvio padrão.
data7193_probes_ensembl <- data_7193[data_7193$rn%in%df_ensembl_duplicated$probes,]
data_7193<-data_7193[!data_7193$rn%in%df_ensembl_duplicated$probes,]

unique_ensemble_duplicated <- unique(df_ensembl_duplicated$ensembl_gene_id)
data7193_sd <- data.frame()
data7193_sd_rbind <- data7193_sd

for(ii in 1:length(unique_ensemble_duplicated)){
  df7193_subset_ensembl <- df_ensembl_duplicated[df_ensembl_duplicated$ensembl_gene_id==unique_ensemble_duplicated[ii],]
  df7193_subset_data <- data7193_probes_ensembl[data7193_probes_ensembl$rn%in%df7193_subset_ensembl$probes,]
  
  df7193_subset_data<- data.frame(df7193_subset_data[,-1],row.names = df7193_subset_data$rn)
  #df_subset_data <- df_subset_data[,-1]
  data7193_sd <- as.data.frame(df7193_subset_data[which.max(rowSds(as.matrix(df7193_subset_data))),])
  data7193_sd_rbind <- rbind(data7193_sd_rbind,data7193_sd)
}

temp7193 <- df_gene_map_aduplicated[isUnique(df_gene_map_aduplicated[,2]),]
temp_data7193 <- data_7193[data_7193$rn%in%temp7193$probes,]

data7193_sd_rbind<-setDT(data7193_sd_rbind,keep.rownames = TRUE)
data7193_rbind_all<-rbind(temp_data7193,data7193_sd_rbind)
data7193_rbind_all<-data.frame(data7193_rbind_all[,-1],row.names = data7193_rbind_all$rn)

#data_7193 <- exprs(exp_7193)

# library(rat2302.db)
# 
# x <- rat2302ENTREZID
# entrez <- unlist(as.list(x[rownames(data)]))
# x <- rat2302SYMBOL
# symbol <- unlist(as.list(x[rownames(data)]))


#colnames(data) <- c(rep("control_10w", 3), rep("diabetic_10w", 3), rep("control_20w", 3), rep("diabetic_20w",3))
tps_7193 <- c(rep("lean_12w", 3), rep("obese_12w", 3))

t_7193 <- factor(tps_7193)

#ph = dt@phenoData
#list = list.files(celpath,full.names=TRUE)

#Dizer a qual grupo determinada amostra pertence
#treats <- as.factor(colnames(data))


design_7193 <- model.matrix(~0+t_7193)

fit_7193 <- lmFit(data7193_rbind_all, design_7193)
#fit$genes$entrez <- entrez
#fit$genes$symbol <- symbol


contrasts_7193 <- makeContrasts(t_7193obese_12w-t_7193lean_12w, levels = design_7193)

ct.fit_7193 <- eBayes(contrasts.fit(fit_7193, contrasts_7193))

res.fit_7193 <- decideTests(ct.fit_7193, method = "global", adjust.method = "BH", p.value = 0.05)

degs_7193 <- data.frame(LogFC = ct.fit_7193$coefficients, p.value = ct.fit_7193$p.value, degenes = unclass(res.fit_7193),stringsAsFactors = F)

ups_7193 <- subset(degs_7193,degs_7193[,3]==1)
downs_7193 <- subset(degs_7193,degs_7193[,3]==-1)
degs_up_down_7193 <- rbind(downs_7193,ups_7193)

probes7193_up_down <- rownames(degs_up_down_7193)
genes7193_up_down <- df_gene_map_aduplicated[df_gene_map_aduplicated$probes%in%probes7193_up_down,]
degs_up_down_7193 <- setDT(degs_up_down_7193,keep.rownames = TRUE)
genes7193_up_down <- merge(degs_up_down_7193,genes7193_up_down,by.x ="rn",by.y = "probes",all.y = TRUE)
colnames(genes7193_up_down)[c(1,2,3,4)] <- c("probes","logFC","pvalue","degs")

degs_universe_7193 <- degs_7193
probes7193_universe <- rownames(degs_universe_7193)
genes7193_universe <- df_gene_map_aduplicated[df_gene_map_aduplicated$probes%in%probes7193_universe,]
degs_universe_7193 <- setDT(degs_universe_7193,keep.rownames = TRUE)
genes7193_universe <- merge(degs_universe_7193,genes7193_universe,by.x ="rn",by.y = "probes",all.y = TRUE)
colnames(genes7193_universe)[c(1,2,3,4)] <- c("probes","logFC","pvalue","degs")

write.csv(genes7193_up_down,file="/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_12wObeseX12wLean.csv",row.names = FALSE)


library(clusterProfiler)
library(rat2302.db)
ego7193 <- enrichGO(gene         = genes7193_up_down$ensembl_gene_id,
                 universe = genes7193_universe$ensembl_gene_id,
                 OrgDb         = rat2302.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)
enrich_data7193 <- ego7193@result
enrich_cutOff7193 <- enrich_data7193[enrich_data7193$p.adjust <= 0.05,]

# Pega a coluna do geneID , trata as barras , da um unique e xablau !

enriched_genes7193 <- enrich_cutOff7193[,8]
enriched_genes7193 <- strsplit(enriched_genes7193,split = "/")
enriched_genes7193 <- unlist(enriched_genes7193)
enriched_genes7193 <- unique(enriched_genes7193)

enriched_genes7193 <- genes7193_up_down[genes7193_up_down$ensembl_gene_id%in%enriched_genes7193]
write.csv(enriched_genes,file="/home/eulle/Documents/Mestrado/MicroArray/results_new_method/genes_12wObeseX12wLean_enriched.csv",row.names = FALSE)

#IDs.up <- ups
#IDs.downs <- rownames(downs)

#library("biomaRt")

# aaa <- topTable(ct.fit, number = 200, adjust.method = "BH")
# 
# library(biomaRt)
# dt_ensembl <- read.table(file = "/home/eulle/Documents/Mestrado/MicroArray/mart_export.txt",header = T,na.strings = "NA" ,sep = ",")
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL")
# mart <- useDataset("rnorvegicus_gene_ensembl",mart)
# #g <- getGene( id = "1370583_s_at",type = "affy_rat230_2", mart = mart)
# annotLookup <- getBM(mart=mart, attributes=c("affy_rat230_2", "ensembl_gene_id", "gene_biotype", "external_gene_name"), filter="affy_rat230_2", values=rownames(data), uniqueRows=TRUE)
# 
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
# genes_ups <- merge(genes_ups,ups,by="row.names",all.x=TRUE)
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
