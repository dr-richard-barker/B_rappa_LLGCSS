.libPaths(c("/tmp/R/library", .libPaths()))
BiocManager::install("pathview")
library(tximport)
library(DESeq2)

library(tidyverse)
library(BiocManager) # Needed to install packages from Bioconductor
library(ggfortify)

library(ComplexHeatmap)

library(EnhancedVolcano)

library(tidyHeatmap)

library(clusterProfiler)

library(pathview)

library(goseq)

library(fgsea)

library(enrichplot)

library(ggnewscale)
dir.create("Brapa_analysis/04-DESeq2_NormCounts/PCA", recursive = TRUE)

dir.create("Brapa_analysis/05-DESeq2_DGE/DGE_plots/PCA", recursive = TRUE)

dir.create("Brapa_analysis/05-DESeq2_DGE/DGE_plots/Heatmaps", recursive = TRUE)

dir.create("Brapa_analysis/05-DESeq2_DGE/DGE_plots/VolcanoPlots", recursive = TRUE)

dir.create("Brapa_analysis/05-DESeq2_DGE/DGE_plots/GSEA", recursive = TRUE)

dir.create("Brapa_analysis/05-DESeq2_DGE/DGE_plots/pathview", recursive = TRUE)
list.dirs(path = "Brapa_analysis", full.names = TRUE, recursive = TRUE) %>% cat(sep = "\n")
counts_dir <- "Brapper_fastq/RSEM_output"

metadata_dir <- "./"

norm_output <- "Brapa_analysis/04-DESeq2_NormCounts"

pca_dir <- "Brapa_analysis/04-DESeq2_NormCounts/PCA"

DGE_output <- "Brapa_analysis/05-DESeq2_DGE"

DGE_pca <- "Brapa_analysis/05-DESeq2_DGE/DGE_plots/PCA"

DGE_heatmaps <- "Brapa_analysis/05-DESeq2_DGE/DGE_plots/Heatmaps"

DGE_volcano <- "Brapa_analysis/05-DESeq2_DGE/DGE_plots/VolcanoPlots"

DGE_gsea <- "Brapa_analysis/05-DESeq2_DGE/DGE_plots/GSEA"

pathview_dir <- "Brapa_analysis/05-DESeq2_DGE/DGE_plots/pathview"



## Get the KEGG annotations

# kegg_organism = "bra"

# brapa_kegg_pathway_names <- read.delim("brapa_kegg_pathway_names.tsv", header=FALSE, col.names=c("pathway_id", "pathway_name"))

# brapa_gene_to_kegg <- read.delim("brapa_gene_to_kegg.tsv", header=FALSE, col.names=c("pathway_id", "gene_id"))

# brapa_symbol_to_kegg_id_map <- read.csv("brapa_symbol_to_kegg_id_map.csv")
sampleTable <- read.csv(file.path(metadata_dir,"Brapa_metadata.csv"),

                  header=TRUE, row.names=1, stringsAsFactors=TRUE, strip.white=TRUE, sep=",")
sampleTable
group <- sampleTable[,1]

group_names <- paste0("(",group,")",sep = "") # define group names

group <- make.names(group) # coerce group names to be compatible with R models

names(group) <- group_names
contrasts <- combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison

contrast.names <- combn(levels(factor(names(group))),2)

contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names

contrasts <- cbind(contrasts,contrasts[c(2,1),])

colnames(contrasts) <- contrast.names
contrasts
## List RSEM raw gene count files ##

files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)



## Reorder the *genes.results files to match the ordering of the samples in the sampleTable ##

files <- files[vapply(rownames(sampleTable), function(x) grep(paste0(counts_dir, '/', x, ".genes.results$"), files, value=FALSE), FUN.VALUE=1)]



names(files) <- rownames(sampleTable)
files
#txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)



## Add 1 to genes with lengths of zero - needed to make DESeqDataSet object ##

#txi.rsem$length[txi.rsem$length == 0] <- 1
#dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~Factor.Value)
#summary(dds)
#keep <- rowSums(counts(dds)) > 10

#dds <- dds[keep,]
#summary(dds)
## Create a data frame containing the RSEM raw counts that were imported with tximport ##

#rawCounts <- as.data.frame(txi.rsem$counts)



## Add 1 to every gene count then log2 transform the raw counts data ##

#exp_raw <- log2(rawCounts+1)



## Calculate the principal components of the raw counts data ##

#PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
#summary(PCA_raw)
#write.csv(PCA_raw$x,file.path(pca_dir,"OSD-104_PCA_raw_table.csv"), row.names=TRUE)
#autoplot(PCA_raw, data=sampleTable, colour='condition',

#         label=TRUE, label.size=5, size=4, alpha=1) + theme_classic(base_size = 16)
#ggsave(file.path(pca_dir,'OSD-104_PCA_raw_wlabels.png'), width = 8.5, height = 6, dpi = 300)
#autoplot(PCA_raw, data=sampleTable, colour='condition',

#         label=FALSE, label.size=5, size=4, alpha=1) + theme_classic(base_size = 16)



## Save your PCA plot without sample labels

#ggsave(file.path(pca_dir,'OSD-104_PCA_raw_nolabels.png'), width = 8.5, height = 6, dpi = 300)
#dds_1 <- estimateSizeFactors(dds)
#dds_1 <- estimateDispersions(dds_1)
#plotDispEsts(dds_1)
#dds_1 <- nbinomWaldTest(dds_1)

#saveRDS(dds_1, file = "Brapa_analysis/dds.rds")

dds_1 <- readRDS("Brapa_analysis/dds.rds")

normCounts = as.data.frame(counts(dds_1, normalized=TRUE))
## Add 1 to every gene count then log2 transform the normalized counts data ##

exp_norm <- log2(normCounts+1)



## Calculate the principal components of the normalized counts data ##

PCA_norm <- prcomp(t(exp_norm), scale = FALSE)
autoplot(PCA_norm, data=sampleTable, colour='condition',

         label=FALSE, label.size=5, size=4, alpha=1) + theme_classic(base_size = 16)



## Save the PCA plot (without sample labels)

ggsave(file.path(pca_dir,'OSD-104_PCA_norm_nolabels.png'), width = 8.5, height = 6, dpi = 300)
DGE_output_table <- normCounts
options(repr.matrix.max.cols=150)

head(DGE_output_table)
contrasts
for (i in 1:dim(contrasts)[2]){

   res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))

	res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]

	colnames(res_1) <-c(paste0("Log2fc_",colnames(contrasts)[i]),

                        paste0("Stat_",colnames(contrasts)[i]),

                        paste0("P.value_",colnames(contrasts)[i]),

                        paste0("Adj.p.value_",colnames(contrasts)[i]))

	DGE_output_table <- cbind(DGE_output_table,res_1)

}
DGE_output_table$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

DGE_output_table$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
tcounts <- as.data.frame(t(normCounts))

tcounts$group <- group

group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))

group_means <- group_means[-c(1),]

colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))



DGE_output_table <- cbind(DGE_output_table,group_means)
group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))

group_stdev <- group_stdev[-c(1),]

colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))



DGE_output_table <- cbind(DGE_output_table,group_stdev)
organism <- "BRASSICA RAPA"

ann.dbi <- "org.Bra.eg.db"

taxonid <- 3711
keytype <- "ENSEMBL"

annot <- data.frame(rownames(DGE_output_table), stringsAsFactors = FALSE)

colnames(annot)[1]<-keytype

if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){

	annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),

                         keys = rownames(DGE_output_table),

                         keytype = keytype,

                         column = "SYMBOL",

                         multiVals = "first")

}

if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){

        annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),

                               keys = rownames(DGE_output_table),

                               keytype = keytype,

                               column = "GENENAME",

                               multiVals = "first")

}

if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){

        annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),

                              keys = rownames(DGE_output_table),

                              keytype = keytype,

                              column = "ENSEMBL",

                              multiVals = "first")

}

if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){

        annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),

                             keys = rownames(DGE_output_table),

                             keytype = keytype,

                             column = "REFSEQ",

                             multiVals = "first")

}

if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){

        annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),

                               keys = rownames(DGE_output_table),

                               keytype = keytype,

                               column = "ENTREZID",

                               multiVals = "first")

}
string_db <- STRINGdb$new(version="12", species=taxonid, score_threshold=0)

string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]

string_map <- string_map[!duplicated(string_map$SYMBOL),]

annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")
head(string_map)
pthOrganisms(PANTHER.db) <- organism

panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")

panther <- na.omit(panther)

annot$GOSLIM_IDS <- panther
options(width = 10000)  # Adjust the width as needed to avoid truncation

head(annot, n=1)
DGE_output_table <- cbind(annot,DGE_output_table)

rownames(DGE_output_table) <- NULL



DGE_output_table$GOSLIM_IDS <- vapply(DGE_output_table$GOSLIM_IDS, paste, collapse = " | ", character(1L))
options(repr.matrix.max.cols=150)

head(DGE_output_table, n=1) # Scroll to the right to see all the columns we've added
write.csv(contrasts,file.path(DGE_output,"contrasts.csv"), row.names=FALSE)

write.csv(DGE_output_table,file.path(DGE_output, "differential_expression.csv"), row.names=FALSE)
exp_norm <- log2(normCounts+1)
dge_p <- DGE_output_table[DGE_output_table$`Adj.p.value_(FLT)v(GC)`<0.05,] # filter by adj p-value

dge_p_flc_up <- dge_p[dge_p$`Log2fc_(FLT)v(GC)`>1,] # filter by lfc > 1

dge_p_flc_down <- dge_p[dge_p$`Log2fc_(FLT)v(GC)`< -1,] # filter by lfc < -1

dge_full <- rbind(dge_p_flc_up, dge_p_flc_down) # merge together



dge_p_flc <- dge_full[rowSums(is.na(dge_full)) != ncol(dge_full), ] # remove rows with NA
exp_dge <- exp_norm[dge_p_flc$ENSEMBL,]

dim(exp_dge)
## Calculate the principal components of the normalized counts data ##

PCA_dge <- prcomp(t(exp_dge), scale = FALSE)



## Create PCA plot without samples labeled

autoplot(PCA_dge, data=sampleTable, colour='condition',

         label=FALSE, label.size=5, size=4, alpha=1) + theme_classic(base_size = 16)



## Save PCA plot without samples labeled

ggsave(file.path(DGE_pca,'OSD-104_PCA_DGE_nolabels.png'), width = 8.5, height = 6, dpi = 300)
exp_dge_scale <- (exp_dge - min(exp_dge)) / (diff(range(exp_dge)))
Heatmap(exp_dge_scale, show_row_names=FALSE, row_title="DEGs (N=772)",

   heatmap_legend_param = list(title="Scaled Expression"))
heatmap = Heatmap(exp_dge_scale, show_row_names=FALSE, row_title="DEGs (N=772)",

   heatmap_legend_param = list(title="Scaled Expression"))
save_pdf(heatmap, file.path(DGE_heatmaps,'OSD-104_heatmap_DGE.pdf'),

        width = 5, height = 5, units = c("in"))
# Volcano plot showing genes differentially expressed in FLT vs GC

EnhancedVolcano(DGE_output_table,

    lab = DGE_output_table$SYMBOL,

    x = 'Log2fc_(scent)v(no_scent)',

    y = 'Adj.p.value_(scent)v(no_scent)',

    title = 'Scent versus No Scent',

    legendLabels=c('NS','|Log2FC| > 1','Adj. p-value < 0.05',

      'Adj. p-value < 0.05 & |Log2FC| > 1'),

    pCutoff = 5e-2,

    FCcutoff = 1,

    pointSize = 3.0,

    labSize = 6.0,

    colAlpha=0.5)



## Save your volcano plot

ggsave(file.path(DGE_volcano,'Brapa_volcano_DGE.png'), width = 6.5, height = 8.5, dpi = 300)



## Create a new column with the KEGG IDs

DGE_output_table$KEGG_ID <- brapa_symbol_to_kegg_id_map$KEGG_ID[match(DGE_output_table$SYMBOL, brapa_symbol_to_kegg_id_map$SYMBOL)]



## Create a named vector of fold changes

foldchanges <- DGE_output_table$Log2fc_(scent)v(no_scent)

names(foldchanges) <- DGE_output_table$KEGG_ID



## Generate the pathview plot

pathview(gene.data=foldchanges, pathway.id="00940", species = kegg_organism, limit = list(gene=max(abs(foldchanges)), cpd=1))
mean_exp_cutoff <- 50

rank_var <- "Stat"

keytype <- "KEGG" ## Change ENSEMBL to TAIR for plant studies ##
## Pull DESeq2 results columns from the DGE table ##

IDs <- DGE_output_table %>% dplyr::select( !!rlang::sym(keytype))

all_mean <- DGE_output_table %>% dplyr::select(All.mean)

log2fc <- DGE_output_table %>% dplyr::select( !!rlang::sym("Log2fc_(FLT)v(GC)"))

stat <- DGE_output_table %>% dplyr::select( !!rlang::sym("Stat_(FLT)v(GC)"))

pvalue <- DGE_output_table %>% dplyr::select( !!rlang::sym("P.value_(FLT)v(GC)"))

padj <- DGE_output_table %>% dplyr::select( !!rlang::sym("Adj.p.value_(FLT)v(GC)"))



## Combine all DESeq2 results columns from the DGE table ##

DGE_res <- cbind(IDs, all_mean, log2fc, stat, pvalue, padj)
DGE_res <- DGE_res %>% dplyr::filter(All.mean > mean_exp_cutoff)
dim(DGE_res)
head(DGE_res)
DGE_res_ranked <- DGE_res %>% dplyr::arrange(desc( !!rlang::sym(paste0(rank_var,"_(scent)v(no_scent)"))))
head(DGE_res_ranked)
gene_list <- DGE_res_ranked %>% dplyr::select( !!rlang::sym(paste0(rank_var,"_(scent)v(no_scent)"))) %>% pull



## Add IDs from the DGE_res_ranked dataframe to the gene_list ##

names(gene_list) <- DGE_res_ranked %>% dplyr::select( !!rlang::sym(keytype)) %>% pull



## View your ranked gene_list ##

head(gene_list)
ont <- "BP"

minGSSize <- 10

maxGSSize <- 500

pvalueCutoff <- 0.05

pAdjustMethod <- "BH"
brapa_kegg_pathway_names <- read.delim("brapa_kegg_pathway_names.tsv", header=FALSE, col.names=c("pathway_id", "pathway_name"))

brapa_gene_to_kegg <- read.delim("brapa_gene_to_kegg.tsv", header=FALSE, col.names=c("pathway_id", "gene_id"))



## Perform GSEA using the KEGG annotations

gse <- GSEA(gene_list, TERM2GENE=brapa_gene_to_kegg, TERM2NAME=brapa_kegg_pathway_names, pvalueCutoff = 0.5)
gse_table <- as.data.frame(gse)



## Save the GSEA output table to a file ##

write.csv(gse_table,file.path(DGE_gsea, paste("Brapa_scent_v_no_scent_GSEA", rank_var, "ranked_output.csv", sep="_")))



head(gse_table)
require(DOSE)

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)



## Save your dotplot ##

ggsave(file.path(DGE_gsea, paste("Brapa_scent_v_no_scent_GSEA", rank_var, "ranked_dotplot.png", sep="_")), width = 11, height = 8.5, dpi = 300)
gse_pairwise_termsim <- pairwise_termsim(gse)

emapplot(gse_pairwise_termsim, showCategory = 10)



## Save your network enrichment map ##

ggsave(file.path(DGE_gsea, paste("Brapa_scent_v_no_scent_GSEA", rank_var, "ranked_network_map.png", sep="_")), width = 11, height = 8.5, dpi = 300)
## print session info ##

print(" ")

print("Session Info below: ")

sessionInfo()
