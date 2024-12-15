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

# Choose the files 
TEcounts = read_table("./data/9_TEtranscripts_unstranded_norep4.25oC/fbxa215_TEtranscripts_withoutrep4_25oC.cntTable")

# manually add names to columns

colnames(TEcounts)[1] <- c("TE") 
  
# keep everything after withoutRep4
colnames(TEcounts)[2:8] <- str_replace(colnames(TEcounts)[2:8],
                                       '.+_withoutRep4//(.+)', '\\1')
colnames(TEcounts)

# to remove everything after replicate number
colnames(TEcounts)[2:8] <- str_replace(colnames(TEcounts)[2:8],
                                       '(.+)Aligned.sortedByCoord.out.bam.+', '\\1')
  
colnames(TEcounts)


###############################################
# make sample info  and TE counts matrix ---- #
###############################################

sampleInfoTE = as.data.frame(colnames(TEcounts[2:8]))
sampleInfoTE$group = sub("_.*", "", sampleInfoTE$`colnames(TEcounts[2:8])`)

colnames(sampleInfoTE) = c("sample", "group")
row.names(sampleInfoTE) = sampleInfoTE$sample

TEcounts_matrix <- as.data.frame(TEcounts)
rownames(TEcounts_matrix) <- TEcounts_matrix$TE
head(TEcounts_matrix)

TEcounts_matrix <- TEcounts_matrix[,2:8]
head(TEcounts_matrix)


#######################################
# Run DESeq2  ---- 
#######################################

ddsObjTE.raw <- DESeqDataSetFromMatrix(countData = TEcounts_matrix,
                                       colData = sampleInfoTE,
                                       design = ~ group)

ddsObjTE <- DESeq(ddsObjTE.raw)


rlogObjTE = rlog(ddsObjTE)

#######################################
# save normalised and not normalised counts  ---- 
#######################################
normTECounts = assay(rlogObjTE)

notNormTECounts = counts(ddsObjTE, normalized = F)


#######################################
# PCA  ---- 
#######################################

plotPCAEx(rlogObjTE,PCx = 1,PCy = 2,cond = "group",ntop = 500, T)
plotPCAEx(rlogObjTE,PCx = 2,PCy = 3,cond = "group",ntop = 500, T)


#######################################
# Scatter plots of replicates and non-replicates   ---- 
#######################################


##  Replicate comparison plots

#scatter plot
pdf("./analysis/figures/scatterPlots_RepsTEs.pdf", height = 15, width = 15)
done = c()
par(mfrow=c(2,2))
for (i in colnames(ddsObjTE)){
  for( j in colnames(ddsObjTE)){
    if( !(paste(i,j,sep="") %in% done)){
      if(i != j){
        if( sub("\\d\\D$","",i) == sub("\\d\\D$","",i)){
          plot_reps(ddsObjTE,x =i ,y = j)
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

### CONTINUE FROM HERE


comp1 = c( "N2")
comp2 = c( "SX3678")




##Now gete the DEgs in a list. FIRST ANALYSIS FOR  0.05 padj
degList_TE = list()
degList2_TE = list()
allDegs_TE = c()

#FDR is alpha on line 226. silencing lines 228-230 to set up a log fold change myself, to not shrink lfc.

for(k in 1:length(comp1)){
  i = comp1[k]
  j = comp2[k]
  name = paste(i,j,sep = "_")
  degList_TE[[name]] = results(ddsObjTE, contrast = c("group", 
                                                     paste(i, sep = "_"), 
                                                     paste(j,sep = "_")),
                                alpha = 0.05)
  #degList_TE[[name]] = lfcShrink(ddsObj, 
  #                        res = degList_TE[[name]],
  #                         type = "ashr")
 
  
  degList2_TE[[name]]= processRes(degList_TE[[name]],FDR = 0.05,lfc = 0)
  summary(degList_TE[[name]])
  allDegs_TE = c(allDegs_TE, row.names(degList2_TE[[name]]))
}

allDegs_TE = unique(allDegs_TE)

length(allDegs_TE)
names(degList2_TE)


# make a counts dataframe only with the TEs. also a degs list only with TEs

normTECounts <-normTECounts[grep(":",rownames(normTECounts)),]
nrow(normTECounts)

length(allDegs_TE)
allDegs_TE <- allDegs_TE[grep(":",allDegs_TE)]
degList_TE[[name]] <- degList_TE[[name]][grep(":",rownames(degList_TE[[name]])),]

# heatmap

myCol = colorRampPalette(brewer.pal(9,"PRGn"))(100)
pheatmap(normTECounts[ allDegs_TE, ],
         col = myCol, scale = "row", show_rownames = T)



#######################################
# MA-plots  ---- 
#######################################



basemeanTE <- degList_TE[[name]]$baseMean
log2foldchange_listTE <- degList_TE[[name]]$log2FoldChange
TE_list <- rownames(degList_TE[[name]])
adjustedpvalue_listTE <- degList_TE[[name]]$padj

length(basemeanTE)
length(log2foldchange_listTE)
length(TE_list)
length(adjustedpvalue_listTE)


forplotsTE <- data.frame(TE_list, basemeanTE, log2foldchange_listTE , adjustedpvalue_listTE)

#View(forplotsTE)

TEs_highlight = as_data_frame(allDegs_TE)
TEs_highlight

forplotsTE_highlight = filter(forplotsTE,
                              TE_list %in% TEs_highlight$value,
                            ignore.case = TRUE)


## write DEseq TE table 

nrow(forplotsTE)
nrow(normTECounts)

normTECounts2 <- as.data.frame(normTECounts)

normTECounts2$TE_list <- row.names(normTECounts2)

row.names(normTECounts2) <- NULL

normTECounts2

TEs_DESeq2_normCounts <- left_join(forplotsTE, normTECounts2, by = "TE_list" )

TEs_DESeq2_normCounts


write.table(TEs_DESeq2_normCounts, file = "TEtranscripts_DESeq2_normCounts_YA_25C.txt", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)


#View(forplotsTE_highlight)

forplotsTE_2 <- anti_join(forplotsTE,
                          forplotsTE_highlight)


forplotsTE_2$highlight = "no"
forplotsTE_highlight$highlight = "yes"
length(forplotsTE_highlight$highlight)

forplotsTE_3 <- bind_rows(forplotsTE_2,forplotsTE_highlight)



maplot_pretty_TE <- ggplot(forplotsTE_3, aes(x = log2(basemeanTE),
                                        y = log2foldchange_listTE, color=highlight)) +
  scale_color_manual(values=c('#000000','#CD2626'), labels = c("> 0.05", "< 0.05")) + 
  geom_point(aes(alpha = highlight)) +
  scale_alpha_manual(values=c(0.3,1)) +
  labs(x = "log2(Mean of Normalized counts)", y = "log2(fold change)",
       title = expression(paste("TE expression in N2 vs ", italic("fbxa-215"))),
       color = "Adj. p-value") +
  theme_bw(base_size = 14)

maplot_pretty_TE + coord_cartesian(ylim = c(-5,5))


### volcano plot

volcanoTEtra <- ggplot(data=forplotsTE_3, aes(x=log2foldchange_listTE, y=-log10(adjustedpvalue_listTE),
                                              color = highlight)) +
  
  geom_hline(yintercept=0, color='#000000', na.rm=TRUE) +
  geom_hline(yintercept=-log10(0.05), color='#000000', na.rm=TRUE, linetype = "dashed") +
  geom_vline(xintercept=0, color='#000000', na.rm=TRUE) +
  geom_vline(xintercept=log2(2), color='#000000', na.rm=TRUE, linetype = "dashed") +
  geom_vline(xintercept=-log2(2), color='#000000', na.rm=TRUE, linetype = "dashed") +
  scale_color_manual(values=c('#000000','#CD2626'), labels = c("> 0.05", "< 0.05")) + 
  scale_fill_manual(values=c('#000000','#CD2626')) + 
  geom_point() +
  labs(x = "log2(fold change)", y = "-log10(p-value)",
       title = expression(paste("TE expression in N2 vs ", italic("fbxa-215"))),
       color = "Adj. p-value") +
  geom_text_repel(
    data = subset(forplotsTE_3, highlight == "yes"),
    aes(label = TE_list),
    show.legend = FALSE,
    size = 3
    #box.padding = unit(0.35, "lines"),
    # point.padding = unit(0.3, "lines")
  ) +
  theme_bw(base_size = 14)

volcanoTEtra + coord_cartesian(xlim = c(-2,2))

































