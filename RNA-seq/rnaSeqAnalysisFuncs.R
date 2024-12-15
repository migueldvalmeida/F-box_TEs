####################################################################################################
# Functions for the analysis of RNASeq data (or data processed with DESeq2)
# source("rnaSeqAnalysisFuncs.R") to get the functions in the environment
# From Jon Price. 
####################################################################################################

####################################################################################################
# REQUIRED LIBRARIES
####################################################################################################
library(genefilter) 
library(ggplot2)
library(lattice)
library(eulerr)

####################################################################################################
#' Scatter log2 (X-Y) plot for quality control  
#'
#' This function takes a condition of choice as the value for cond_choice
#' and the replicates of choice x and y 
#' This function can be embedded in a for loop to plot all samples against
#' eachother. 
#' 
#' @param x The first column of the DDS matrix to be plotted
#' @param y The second column of the DDS matrix to be plotted
#' @return No return but plots the x-y scatterplot
#' @examples
#' plot_reps(dds, 1,2,1,"group");
#' @export

plot_reps =  function(dds,x=1,y=2){
  ##  Estimate the size factors for normalisation
  dds<-estimateSizeFactors(dds)
  
  ## Extract the normalised counts for the condition you want
  rep_values<- counts(dds, normalized=TRUE)[,c(x,y)]
  
  # Take logs of these values
  vals <- log2(rep_values[,c(x,y)] + 1)
  # And plot
  plot(vals,pch=16, cex=0.4,xlab=paste('rep',x),ylab=paste('rep',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  #title(paste("Comparison of",cond_choice,"replicates"))
}



####################################################################################################
#' Modified plotPCA 
#  Add ability to select principal components to be plotted (PCx/PCy)
#  Optional labelling of points, based on http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html
#' 
#' @param object This object has to be a tranformed dds object (rlog() or vst() )
#' @param PCx The First principle component to be plotted
#' @param PCy The second principle component to be plotted 
#' @param cond To colour the plots choose a column from colData(object) (this can be a vector for combinatorial use)
#' @param ntop The cutoff - ntop variance genes are selected 
#' @param labels T or F plot sample labels or not
#' @return No return but plots the PCA plot
#' @examples
#' plotPCAEx(vstDDs, 1, 2, "group", 1000, T)
#' @export

plotPCAEx = function(object, PCx = 1, PCy = 2, cond="condition", ntop=500, labels = TRUE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(cond %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  cond.df <- as.data.frame(colData(object)[, cond, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(cond) > 1) {
    factor(apply( cond.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[cond]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PCa=pca$x[,PCx], PCb=pca$x[,PCy], group=group, cond.df, name=colnames(object))
  
  pc1 <- ggplot(data=d, aes_string(x="PCa", y="PCb", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC",PCx,": ",round(percentVar[PCx] * 100),"% variance")) +
    ylab(paste0("PC",PCy,": ",round(percentVar[PCy] * 100),"% variance")) +
    coord_fixed()
  
  pc2 <- pc1 + geom_point()
  pc3 <- pc2 +     theme_classic() 

  #  Finally add the labels, using ggrepel so that they dont write over each other or the points  
  if (labels)
  {
    library("ggrepel")
    pc3 + geom_text_repel(aes(label = name),
                    color = "gray20",
                    data = d,
                    force = 10)
  } else {pc3}
}



####################################################################################################
#' Process results object created by DESeq2 "results()" 
#' Commonly you want to filter the results object to only contain the "significant" genes
#' This function returns a subsetted data-frame containing only the significant genes
#' 
#' @param red This object has to be a tranformed dds object (rlog() or vst() )
#' @param FDR The FDR cutoff
#' @param lfc The logFoldChange cutoff 
#' @return returns a subsetted Res object
#' @examples
#' subsetRes = processRes(degs, 0.01, 1)
#' @export

processRes = function(res, FDR, lfc){
  print(summary(res))
  resOriginal = res
  res = res[order(res$pvalue),]
  #get table of significant DEGs
  res = subset(res, !is.na(padj))
  res = res[res$padj < FDR & sqrt((res$log2FoldChange)^2) > lfc,]
  print(summary(res))
  return(res)
}



####################################################################################################
#' Plots a 3 way Euler diagram based on 3 gene sets 
#' @param anames Genes for set a (vector)
#' @param bnames Genes for set b (vector)
#' @param cnames Genes for set c (vector)
#' @param a Sample name for set a (string)
#' @param b Sample name for set b (string)
#' @param c Sample name for set c (string)
#' @return returns list of length 2 - plot() of the 1st will plot the diagram and the gene lists are contained in 2
#' @examples
#' eulerList = plotEuler(degsSample1,degsSample2,degsSample3, "sample1", "sample2", "sample3")
#' @export

plotEuler = function(anames,bnames,cnames,a,b,c){
  ab = paste(a,"&",b, sep = "") 
  ac = paste(a,"&",c, sep = "") 
  bc = paste(b,"&",c, sep = "") 
  abc= paste(a,"&",b,"&",c, sep = "") 
  
  lena = length(anames)
  lenb = length(bnames)
  lenc = length(cnames)
  lenab = length(which(anames %in% bnames))
  lenac = length(which(anames %in% cnames))
  lenbc = length(which(bnames %in% cnames))
  lenall =  length(which(anames[which(anames %in% bnames)] %in% anames[which(anames %in% cnames)]))
  
  geneList = list() 
  geneList[[abc]]   = anames[which(anames[(anames %in% bnames)] %in% cnames)]
  geneList[[a]]   = anames[!(anames %in% bnames) & !(anames %in% cnames)] 
  geneList[[b]]   = bnames[!(bnames %in% anames) & !(bnames %in% cnames)] 
  geneList[[c]]   = cnames[!(cnames %in% bnames) & !(cnames %in% anames)] 
  geneList[[ab]]   = anames[(anames %in% bnames)][!( anames[(anames %in% bnames)] %in%  cnames) ]
  geneList[[ac]]   = anames[(anames %in% cnames)][!( anames[(anames %in% cnames)] %in%  bnames) ]
  geneList[[bc]]   = cnames[(cnames %in% bnames)][!( cnames[(cnames %in% bnames)] %in%  anames) ]
  
  
  
  
  fit1 <- euler( structure(c(lena - lenab - lenac + lenall,
                             lenb - lenab - lenbc + lenall,
                             lenc - lenac - lenbc + lenall,
                             lenab - lenall,
                             lenac - lenall,
                             lenbc - lenall,
                             lenall), names = c(a,b,c,ab,ac,bc,abc)))
  
  #check residuals
  dotplot(resid(fit1), xlab = "",
          panel = function(...) {
            panel.abline(v = 0, lty = 2)
            panel.dotplot(...)
          })
  error_plot(fit1)
  
  #Just check the coefficiants
  coef(fit1)
  return(list(fit1, geneList)) 
}







