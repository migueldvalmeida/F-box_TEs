## libraries needed:
library(RColorBrewer)
library(Gviz)
library(GenomicFeatures)
library(data.table)
library(tidyverse)

gen <- 'WBcel235'

### TWEAKED BECAUSE THE TRACKS FROM THE PAPER ARE IN UCSC NAMES
## gene coordinates
mygenes <- fread("genes_to_plot.txt")
gene <- mygenes[7] ## change according to which gene on the table you want to plot

chr <- gene[["chromosome"]]
mygene_start <- gene[["start"]]
mygene_end <- gene[["end"]]
mygene <- gene[["name"]]

# create track with gene annotations:
txdb <- makeTxDbFromGFF("A2_fbx.gtf", format="gtf")
options(ucscChromosomeNames=FALSE)
grtrack_A2 <- GeneRegionTrack(
  txdb,
  genome=gen,
  chromosome=chr,
  name='A2',
  collapseTranscripts=FALSE,
  shape="box",
  stacking="dense",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#e52521",
  col= NULL,
  col.line=NULL
)


txdb2 <- makeTxDbFromGFF("A1_fbx.gtf", format="gtf")
options(ucscChromosomeNames=FALSE)
grtrack_A1 <- GeneRegionTrack(
  txdb2,
  genome=gen,
  chromosome=chr,
  name='A1',
  collapseTranscripts=FALSE,
  shape="box",
  stacking="dense",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#010101",
  col= NULL,
  col.line=NULL
)


txdb3 <- makeTxDbFromGFF("B_fbx.gtf", format="gtf")
options(ucscChromosomeNames=FALSE)
grtrack_B <- GeneRegionTrack(
  txdb3,
  genome=gen,
  chromosome=chr,
  name='B',
  collapseTranscripts=FALSE,
  shape="box",
  stacking="dense",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#2b4b9b",
  col= NULL,
  col.line=NULL
)


## adapted the ce11_rmsk_TE.gtf file to non-UCSC chr names by doing: sed 's/^chr//' ce11_rmsk_TE.gtf > ce11_rmsk_TE_forTracks.gtf in command line

txdb4 <- makeTxDbFromGFF("ce11_rmsk_TE_forTracks.gtf", format="gtf")
options(ucscChromosomeNames=FALSE)
TEtrack <- GeneRegionTrack(
  txdb4,
  genome=gen,
  chromosome=chr,
  name='TEs',
  collapseTranscripts=FALSE,
  shape="box",
  stacking="dense",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#000000",
  col= NULL,
  col.line=NULL
)



axisTrack <- GenomeAxisTrack()


track_list <- list(axisTrack,
                   grtrack_A1,
                   grtrack_A2,
                   grtrack_B,
                   TEtrack)

plotTracks(track_list,
           chromosome=chr,
           from = mygene_start + 400000,
           to = mygene_end - 200000,
           background.title="white",
           col.title="black",
           col.axis="black")


track_list2 <- list(axisTrack,
                    grtrack_A1,
                    grtrack_A2,
                    grtrack_B)

plotTracks(track_list2,
           chromosome=chr,
           from = mygene_start + 400000,
           to = mygene_end - 200000,
           background.title="white",
           col.title="black",
           col.axis="black")




