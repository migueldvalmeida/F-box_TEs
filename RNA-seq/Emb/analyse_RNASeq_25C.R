#######################################
#Load libraries  ---- 
#######################################

library(ggplot2)
library(pheatmap)
library(reshape2)
library(genefilter) 
library(ggrepel)
library(eulerr)
library(DESeq2)
library(UpSetR)
library(GenomicFeatures)
library(biomaRt)
library(tidyverse)
library(tximport)
library(RColorBrewer)
library(ashr)
library(ggrepel)


#######################################
#Load Counts  ---- 
#######################################

# Get the tx2gene table
tx2gene = read.table("./info/tx2gene.tsv", header = T)

# Choose the files 
files = list.files(pattern = ".*sf", 
                   path =  "./4_salmon_count/", 
                   recursive = T,
                   full.names = T)
# add names
names(files) = list.files(path =  "./4_salmon_count/")


# Remove a sample if nessecary 
#files = files[-12]

# scaled using the average transcript length over samples and then the library size (lengthScaledTPM)
#nice 
tpm <- tximport(files, type = "salmon",
                tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")


# This is the equivelant of raw counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "no")


#######################################
# make sample info  ---- 
#######################################

sampleInfo = as.data.frame(colnames(txi$counts))
sampleInfo$group = sub("_.*", "", sampleInfo$`colnames(txi$counts)`)

colnames(sampleInfo) = c("sample", "group")
row.names(sampleInfo) = sampleInfo$sample



#######################################
# Make annotation   ---- 
#######################################

# See what annotations are available
listMarts()

# Set ENSEMBL as the annotation and see what species are there
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets, 50)

# choose c.elegans
ensembl = useDataset("celegans_gene_ensembl",mart=ensembl)

# Filters are the info we have 
filters = listFilters(ensembl)
filters[1:5,]

# attributes are the information we want to get 
attributes = listAttributes(ensembl)
attributes[1:50,]


# Get the annotation 
annotation = getBM(attributes=c('ensembl_gene_id',
                                'description' , 'gene_biotype', 
                                'external_gene_name','external_gene_source',
                                'entrezgene_accession',"chromosome_name","start_position"), 
                   filters = 'ensembl_gene_id', 
                   values = row.names(txi$abundance), 
                   mart = ensembl, uniqueRows = T)

# See if any have duplicate entries 
annotation[duplicated(annotation$ensembl_gene_id),]
annotation
dim(annotation)


# isolate protein coding genes
protein_codingGenes = annotation[annotation$gene_biotype == "protein_coding",]
nrow(protein_codingGenes)


# If you want to add other annotation sets:
#t1 = read.table("../info/morc1/exampleGeneSet.txt", header = F, sep = "\t")
#protein_codingGenes$example = 0
#protein_codingGenes[protein_codingGenes$ensembl_gene_id %in% t1$V1,]$example = 1




#######################################
# Run DESeq2  ---- 
#######################################

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleInfo,
                                       design = ~ group)
ddsObj.raw = ddsObj.raw[protein_codingGenes$ensembl_gene_id,]

ddsObj <- DESeq(ddsObj.raw)

plotDispEsts(ddsObj)

rlogObj = rlog(ddsObj)

#######################################
# save normalised and not normalised counts  ---- 
#######################################
normCounts = assay(rlogObj)
normCountsAnnot = cbind(protein_codingGenes, normCounts )

notNormCounts = counts(ddsObj, normalized = F)
notNormCountsAnnot = cbind(protein_codingGenes, notNormCounts  )


#######################################
# PCA  ---- 
#######################################

plotPCAEx(rlogObj,PCx = 1,PCy = 2,cond = "group",ntop = 500, T)
plotPCAEx(rlogObj,PCx = 2,PCy = 3,cond = "group",ntop = 500, T)
plotPCAEx(rlogObj,PCx = 1,PCy = 3,cond = "group",ntop = 500, T)


#######################################
# Scatter plots of replicates and non-replicates   ---- 
#######################################


##  Replicate comparison plots

#scatter plot
pdf("./analysis/figures/scatterPlots_Reps.pdf", height = 15, width = 15)
done = c()
par(mfrow=c(2,2))
for (i in colnames(ddsObj)){
    for( j in colnames(ddsObj)){
        if( !(paste(i,j,sep="") %in% done)){
            if(i != j){
                if( sub("\\d\\D$","",i) == sub("\\d\\D$","",i)){
                    plot_reps(ddsObj,x =i ,y = j)
                    abline(a = 0,b = 1)
                    done = c(done, paste(i,j,sep=""),paste(j,i,sep=""))
                    #print(done)
                }
            }
        }
    }
}
dev.off()





#######################################
# DEGs  ---- 
#######################################



comp1 = c( "N2")
comp2 = c( "SX3678")




##Now gete the DEgs in a list. FIRST ANALYSIS FOR  0.05 padj
degList_005 = list()
degList2_005 = list()
allDegs_005 = c()

#FDR is alpha on line 226. silencing lines 228-230 to set up a log fold change myself, to not shrink lfc.

for(k in 1:length(comp1)){
    i = comp1[k]
    j = comp2[k]
    name = paste(i,j,sep = "_")
    degList_005[[name]] = results(ddsObj, contrast = c("group", 
                                                   paste(i, sep = "_"), 
                                                   paste(j,sep = "_")),
                              alpha = 0.05)
    #degList_005[[name]] = lfcShrink(ddsObj, 
        #                        res = degList_005[[name]],
       #                         type = "ashr")
    degList_005[[name]]$geneName = protein_codingGenes$external_gene_name
    degList_005[[name]]$geneName = protein_codingGenes$chromosome_nam
    degList_005[[name]]$geneName = protein_codingGenes$start_position
    
    degList2_005[[name]]= processRes(degList_005[[name]],FDR = 0.05,lfc = 1)
    summary(degList_005[[name]])
    allDegs_005 = c(allDegs_005, row.names(degList2_005[[name]]))
}

allDegs_005 = unique(allDegs_005)

length(allDegs_005)
names(degList2_005)

myCol = colorRampPalette(brewer.pal(9,"PRGn"))(100)
pheatmap(normCounts[ allDegs_005, ],
         col = myCol, scale = "row", show_rownames = F)




#######################################
# MA-plots  ---- 
#######################################

#View(normCountsAnnot)
#View(degList_005[[name]])

basemean <- degList_005[[name]]$baseMean
log2foldchange_list <- degList_005[[name]]$log2FoldChange
genenames_list <- rownames(degList_005[[name]])
adjustedpvalue_list <- degList_005[[name]]$padj
wormgenenames <- normCountsAnnot[,4]

length(basemean)
length(log2foldchange_list)
length(genenames_list)
length(adjustedpvalue_list)
length(wormgenenames)


forplots <- data.frame(genenames_list, basemean, log2foldchange_list , adjustedpvalue_list,wormgenenames)

#View(forplots)

genes_highlight = as_data_frame(allDegs_005)
genes_highlight

forplots_highlight = filter(forplots,
                               genenames_list %in% genes_highlight$value,
                               ignore.case = TRUE)

#View(forplots_highlight)

forplots_2 <- anti_join(forplots,
                        forplots_highlight)


forplots_2$highlight = "no"
forplots_highlight$highlight = "yes"
length(forplots_highlight$highlight)

forplots_3 <- bind_rows(forplots_2,forplots_highlight)



#remove NAs
length(forplots_3$basemean)
forplots_4 <- na.omit(forplots_3)
length(forplots_4$basemean)

maplot2 <- ggplot(forplots_4, aes(x = log2(basemean),
                                 y = log2foldchange_list, color=highlight)) +
    scale_color_manual(values=c('#000000','#CD2626')) + 
    geom_point(aes(alpha = highlight)) +
    scale_alpha_manual(values=c(0.3,1)) +
    labs(x = "log2(Mean of Normalized counts)", y = "log2(fold change)") +
    geom_text_repel(
        data = subset(forplots_4, highlight == "yes"),
        aes(label = wormgenenames),
        show.legend = FALSE,
        size = 3,
        fontface = "bold") +
    theme_bw()

maplot2 + coord_cartesian(ylim = c(-14,14))







