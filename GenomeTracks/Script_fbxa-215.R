## libraries needed:
library(RColorBrewer)
library(Gviz)
library(GenomicFeatures)
library(data.table)

getwd()

gen <- 'WBcel235'


## gene coordinates
mygenes <- fread("genes_to_plot.txt")
gene <- mygenes[1] ## change according to which gene on the table you want to plot

chr <- gene[["chromosome"]]
mygene_start <- gene[["start"]]
mygene_end <- gene[["end"]]
mygene <- gene[["name"]]

# create track with gene annotations:
txdb <- makeTxDbFromGFF("LongGenes_plus_Transposons_WBcel235.gtf", format="gtf")
options(ucscChromosomeNames=FALSE)
grtrack <- GeneRegionTrack(
  txdb,
  genome=gen,
  chromosome=chr,
  name='WBcel235 genes',
  collapseTranscripts="longest",
  shape="smallArrow",
  stacking="squish",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#000000",
  col= NULL,
  col.line=NULL
)


bw1 <- "N2_merged_reps.bw"
bw2 <- "SX3678_merged_reps.bw"



name1 <- "CPM in N2"
name2 <- "CPM in fbxa-215"




bw1Track <- DataTrack(
  range = bw1,
  genome = gen,
  chromosome = chr,
  name = name1,
  type = "histogram",
  col.histogram="#3f007d",
  fill="#3f007d")


bw2Track <- DataTrack(
  range = bw2,
  genome = gen,
  chromosome = chr,
  name = name2,
  type = "histogram",
  col.histogram="#9e9ac8",
  fill="#9e9ac8")





## adds position in chromosome
axisTrack <- GenomeAxisTrack()


track_list <- list(
  bw1Track,
  bw2Track,
  axisTrack,
  grtrack
)

plotTracks(track_list,
           chromosome=chr,
           from = mygene_start - 500,
           to = mygene_end + 500,
           background.title="white",
           col.title="black",
           col.axis="black",
           ylim=c(0,40))

