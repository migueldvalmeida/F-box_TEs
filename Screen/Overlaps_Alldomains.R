library(tidyverse)
library(taxize)
library(taxizedb)
library(tidyr)
library(taxonomizr)
library(myTAI)
library(dplyr)
library(ggplot2)
library(wordcloud2)
library(RColorBrewer)

###############################################################################
########## ---- Load tables.    ###############################################
###############################################################################

AllDoms <- read.table(file = "Overlaps_allDoms_together.txt", sep = "\t",
                          header = FALSE, fill =TRUE, quote = NULL)

View(AllDoms)
nrow(AllDoms)

# One line got messed up, the line with "M7BBS3" and the next. fix them manually

AllDoms[grepl("unreviewed",AllDoms$V10),] # just one line is the problem

AllDoms[46582,9] = "IPR017981"
AllDoms[46582,9]

# fixed entry of next line for: A0A0F7ZW32

AllDoms[46583,] = c("A0A0F7ZW32","unreviewed","Reverse transcriptase domain-containing protein","1043627","Hirsutella minnesotensis 3608","317","IPR000477","1..154","IPR029021","")
AllDoms[46583,]

# remove last column
AllDoms <- AllDoms[,1:9]

length(unique(AllDoms$V1))
# 40827 unique hits

nrow(AllDoms)

# there's a total of 40827 unique proteins (out of 56559) with TE/viral domains fused to ubiquitous domains.


###############################################################################
########## ---- Add additional taxonomic information to the dataframe.    #####
###############################################################################


# prepared database of taxonomizr, only need to do it once: prepareDatabase('accessionTaxa.sql')


# Extract taxonomy information from classification results
Taxonomy_AllDomsV4 <-getTaxonomy(AllDoms$V4,'accessionTaxa.sql')

#View(Taxonomy_AllDomsV4)

nrow(Taxonomy_AllDomsV4)
nrow(AllDoms)

AllDoms_withTax <- cbind(AllDoms,Taxonomy_AllDomsV4)
#View(AllDoms_withTax)
nrow(AllDoms_withTax)

is.data.frame(AllDoms_withTax)


###############################################################################
########## ---- Remove metagenome and prokaryotes    ##########################
###############################################################################

# remove prokaryotes and viruses 
unique(AllDoms_withTax$superkingdom)

nrow(AllDoms_withTax)
AllDoms_withTaxFiltered <- AllDoms_withTax[!grepl("Bacteria",AllDoms_withTax$superkingdom),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$superkingdom)

AllDoms_withTaxFiltered <- AllDoms_withTaxFiltered[!grepl("Archaea",AllDoms_withTaxFiltered$superkingdom),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$superkingdom)

AllDoms_withTaxFiltered <- AllDoms_withTaxFiltered[!grepl("Viruses",AllDoms_withTaxFiltered$superkingdom),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$superkingdom)

#View(AllDoms_withTaxFiltered)

# What are the superkingdom NAs?
AllDoms_withTaxFiltered_NAs <- AllDoms_withTaxFiltered[is.na(AllDoms_withTaxFiltered$superkingdom),]
#View(AllDoms_withTaxFiltered_NAs)

# most of these are metagenomes... remove
nrow(AllDoms_withTaxFiltered_NAs) # 271 lines with metagenomes and uncultured organisms

AllDoms_withTaxFiltered <- AllDoms_withTaxFiltered[!grepl("metagenome",AllDoms_withTaxFiltered$V5),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$superkingdom)

# removed 237 metagenomes

AllDoms_withTaxFiltered <- AllDoms_withTaxFiltered[!grepl("uncultured",AllDoms_withTaxFiltered$V5),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$superkingdom)

# removed 4 uncultured organisms. 

AllDoms_withTaxFiltered_NAs2 <- AllDoms_withTaxFiltered[is.na(AllDoms_withTaxFiltered$superkingdom),]
#View(AllDoms_withTaxFiltered_NAs2)

# the rest are taxonomy entries with NAs in all fields: prevotella,Stenotrophomonas,Clostridium roseum,
# Endozoicomonas, Pseudonocardia, Lysobacter, Weizmannia, Bradyrhizobium, (bacteria) and Lachancea, and Fomitopsis (fungi)
unique(AllDoms_withTaxFiltered_NAs2$V5)
nrow(AllDoms_withTaxFiltered_NAs2)

# filter out all the NAs, lose only 4 entries of these weird fungi without taxonomy classification... better for simplicity, so I dont have to input data manually.

AllDoms_withTaxFiltered <- AllDoms_withTaxFiltered[!is.na(AllDoms_withTaxFiltered$superkingdom),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$superkingdom)

# this removed 30 NA rows corresponding to the organisms mentioned in lines 111 and 112.  


## need to make a column for Kingdom,so manually add the correspondence between kingdom and phylums.
# removed all the unicellular eukaryotes and kept only multicelullar eukaryotes.

unique(AllDoms_withTaxFiltered$phylum)
write.table(unique(AllDoms_withTaxFiltered$phylum), file = "kingdomToPhylum_original.txt" , col.names = FALSE, row.names = FALSE, quote = FALSE)
# then added the equivalent kingdoms and loaded table again. 


KingdomToPhylum <- read.table(file = "kingdomToPhylum.txt", sep = "\t",
                      header = TRUE, quote = NULL)

#View(KingdomToPhylum)

AllDoms_withTaxFiltered <- left_join(AllDoms_withTaxFiltered,KingdomToPhylum, by = c("phylum" = "Phylum"))
#View(AllDoms_withTaxFiltered)

## then for simplicity also removed all hits in unicellular eukaryotes (2334 hits removed):
unique(AllDoms_withTaxFiltered$Kingdom)
nrow(AllDoms_withTaxFiltered)

AllDoms_withTaxFiltered <- AllDoms_withTaxFiltered[!grepl("UnicelEuk",AllDoms_withTaxFiltered$Kingdom),]
nrow(AllDoms_withTaxFiltered)
unique(AllDoms_withTaxFiltered$Kingdom)



# Add interpro domain description
interpro_TEViral <- read.table("TE_IPR_tofilter_MVA_final.txt", header = TRUE, sep="\t", fill = TRUE, quote = NULL)
interpro_TEViral <- interpro_TEViral[,1:2]
#View(interpro_TEViral)


AllDoms_withTaxFiltered_final <- left_join(AllDoms_withTaxFiltered,interpro_TEViral, by = c("V7" = "Accession"))
#View(AllDoms_withTaxFiltered_final)

# remove quotes from description
 
AllDoms_withTaxFiltered_final$Name <- gsub("\"","",AllDoms_withTaxFiltered_final$Name)
nrow(AllDoms_withTaxFiltered_final)

## add the interpro domain description of the ubiquitous domain. 


ubiquitousDomainList <- read.table("ListOfLargeProteinDomainFamilies_CURATED.txt", header = TRUE, sep="\t", fill = TRUE, quote = NULL)
ubiquitousDomainList <- ubiquitousDomainList[,1:2]
#View(ubiquitousDomainList)



AllDoms_withTaxFiltered_final <- left_join(AllDoms_withTaxFiltered_final,ubiquitousDomainList, by = c("V9" = "Interpro"))
View(AllDoms_withTaxFiltered_final)

# remove quotes from description

AllDoms_withTaxFiltered_final$Domain.description <- gsub("\"","",AllDoms_withTaxFiltered_final$Domain.description)
nrow(AllDoms_withTaxFiltered_final)


# finally edit column names:

AllDoms_withTaxFiltered_final
colnames(AllDoms_withTaxFiltered_final) <- c("UniprotID","status","proteinDescription","taxID","Species","proteinLength",
                                             "TEViralDomainID","TEdomain_location","ubiquitousDomainID","superkingdom",
                                             "phylum","class","order","family","genus","species","Kingdom","TEviralDomain","UbiquitousDomain")
AllDoms_withTaxFiltered_final

## reorder into something more comprehensible and remove the redundant species entry and the status entry

AllDoms_withTaxFiltered_final <- AllDoms_withTaxFiltered_final[,c(1,3,19,9,18,7,8,6,4,10,17,11:16)]

# found two weird entries with not protein domain... weird, removed these manual.

nrow(AllDoms_withTaxFiltered_final)
AllDoms_withTaxFiltered_final <- AllDoms_withTaxFiltered_final[!is.na(AllDoms_withTaxFiltered_final$TEviralDomain),]
nrow(AllDoms_withTaxFiltered_final)

## 26540 hits in multicellular eukaryotes! ----

nrow(AllDoms_withTaxFiltered_final)
length(unique(AllDoms_withTaxFiltered_final$UniprotID))

## 17342 unique proteins, and 9198 that have more than 1 domain in the list ----


# any other NAs in the dataframe? 


na_rows <- AllDoms_withTaxFiltered_final[!complete.cases(AllDoms_withTaxFiltered_final), ]
#View(na_rows)
# all good, all the remaining NAs are from taxonomny columns, from classification ambiguities. e.g. for rotiferans and a fish. 

write.table(AllDoms_withTaxFiltered_final, file = "AllDoms_withTaxFiltered_final.txt" , col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")




#### now manually removing the RNase H domains hits from the RNaseH domains which are too broad and including not only TE-associated proteins.

AllDoms_withTaxFiltered_final_noRNaseH <- AllDoms_withTaxFiltered_final
nrow(AllDoms_withTaxFiltered_final_noRNaseH)

AllDoms_withTaxFiltered_final_noRNaseH <- AllDoms_withTaxFiltered_final_noRNaseH[!AllDoms_withTaxFiltered_final_noRNaseH$TEViralDomainID == "IPR012337",]
nrow(AllDoms_withTaxFiltered_final_noRNaseH) # removed 6651 proteins

AllDoms_withTaxFiltered_final_noRNaseH <- AllDoms_withTaxFiltered_final_noRNaseH[!AllDoms_withTaxFiltered_final_noRNaseH$TEViralDomainID == "IPR022892",]
nrow(AllDoms_withTaxFiltered_final_noRNaseH) # removed 0 proteins

AllDoms_withTaxFiltered_final_noRNaseH <- AllDoms_withTaxFiltered_final_noRNaseH[!AllDoms_withTaxFiltered_final_noRNaseH$TEViralDomainID == "IPR002156",]
nrow(AllDoms_withTaxFiltered_final_noRNaseH) # removed 995 proteins

# final: 18894 protein hits

nrow(AllDoms_withTaxFiltered_final_noRNaseH)
length(unique(AllDoms_withTaxFiltered_final_noRNaseH$UniprotID))

# 12803 unique proteins


write.table(AllDoms_withTaxFiltered_final_noRNaseH, file = "AllDoms_withTaxFiltered_final_noRnaseH.txt" , col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# hits per domain

countUbiqDomain <- AllDoms_withTaxFiltered_final_noRNaseH %>% count(ubiquitousDomainID)
countUbiqDomain

###############################################################################
########## ---- Plotting   ####################################################
###############################################################################

# plot number of single proteins found in each group.

AllDoms_withTaxFiltered_final_noRNaseH_singleProteins <- AllDoms_withTaxFiltered_final_noRNaseH
AllDoms_withTaxFiltered_final_noRNaseH_singleProteins <- AllDoms_withTaxFiltered_final_noRNaseH_singleProteins[!(duplicated(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins[,"UniprotID"])),]
nrow(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins)
nrow(AllDoms_withTaxFiltered_final_noRNaseH)

AllDoms_withTaxFiltered_final_noRNaseH_singleProteins$Kingdom <- ordered(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins$Kingdom,
                                                                 levels = c("Viridiplantae","Metazoa","Fungi"))

plotSingleProteins <- ggplot(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins, aes(x = Kingdom, width=0.70, fill = Kingdom)) +
  geom_bar(aes(y = (..count..))) + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='count', aes(label=..count..), hjust=-0.5) +
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#54278f")) + 
  ylim(0,8500) +
  coord_flip() +
  theme_bw()
plotSingleProteins


# 8.71 x 6.42 

# plot number of domains found in each group.

AllDoms_withTaxFiltered_final_noRNaseH$Kingdom <- ordered(AllDoms_withTaxFiltered_final_noRNaseH$Kingdom,
                                                                         levels = c("Viridiplantae","Metazoa","Fungi"))

plotNumberDomainsCaptured<- ggplot(AllDoms_withTaxFiltered_final_noRNaseH, aes(x = Kingdom, width=0.70, fill = Kingdom)) +
  geom_bar(aes(y = (..count..))) + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='count', aes(label=..count..), hjust=-0.5) +
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#54278f")) + 
  ylim(0,12000) +
  coord_flip() +
  theme_bw()
plotNumberDomainsCaptured

# 8.71 x 6.42 


### to get the word clouds, chose not to do it with single proteins, because then the duplicates have different co-occurring TE domains which wouldnt be represented. 


fungi <- filter(AllDoms_withTaxFiltered_final_noRNaseH, Kingdom == "Fungi")

fungi2 <- fungi %>% count(TEviralDomain)
#View(fungi2)

wordcloud2(fungi2, size = 0.2)

# saved as 836 x 616


top_proteins_fungi2 <- fungi2 %>%
  arrange(desc(n)) %>%  
  slice(1:10)           

top_proteins_fungi2

wordcloud2(top_proteins_fungi2, size = 0.2,ellipticity = 0.3)




metazoa <- filter(AllDoms_withTaxFiltered_final_noRNaseH, Kingdom == "Metazoa")

metazoa2 <- metazoa %>% count(TEviralDomain)
#View(metazoa2)

wordcloud2(metazoa2, size = 0.2)

# saved as 836 x 616

top_proteins_metazoa2 <- metazoa2 %>%
  arrange(desc(n)) %>%  
  slice(1:10)           


top_proteins_metazoa2

wordcloud2(top_proteins_metazoa2, size = 0.2,ellipticity = 0.3)




viridiplantae <- filter(AllDoms_withTaxFiltered_final_noRNaseH, Kingdom == "Viridiplantae")

viridiplantae2 <- viridiplantae %>% count(TEviralDomain)
#View(viridiplantae2)

wordcloud2(viridiplantae2, size = 0.3)

# saved as 836 x 616


top_proteins_viridiplantae2 <- viridiplantae2 %>%
  arrange(desc(n)) %>%  
  slice(1:10)           

top_proteins_viridiplantae2

wordcloud2(top_proteins_viridiplantae2, size = 0.1, ellipticity = 0.3)




### Count number of TE domains per protein 


AllDoms_withTaxFiltered_final_noRNaseH

freq_proteinDom <- table(AllDoms_withTaxFiltered_final_noRNaseH$UniprotID)
#View(freq_proteinDom)

summary_df <- as.data.frame(freq_proteinDom)
colnames(summary_df) <- c("ProteinID", "Frequency")

summary_df <- summary_df[order(summary_df$Frequency, decreasing = TRUE), ]
print(summary_df)

ggplot(summary_df, aes(x = Frequency, y = ..count..)) +
  geom_bar() +
  labs(x = "Number of TE/viral domains", y = "Number of Proteins") +
  ggtitle("Number of TE/viral domains per protein") +
  theme_bw()


# 5.20 x 4.97




###############################################################################
########## ---- Plotting most co-occurring domains  ###########################
###############################################################################


AllDoms_withTaxFiltered_final_noRNaseH

coOccuringDomains <- AllDoms_withTaxFiltered_final_noRNaseH %>% count(ubiquitousDomainID,TEViralDomainID)
coOccuringDomains
#View(coOccuringDomains)

coOccuringDomains$domains <- paste(coOccuringDomains$ubiquitousDomainID, coOccuringDomains$TEViralDomainID, sep = "-")
coOccuringDomains
#View(coOccuringDomains)

# order based on ranking

coOccuringDomains <- coOccuringDomains[order(-coOccuringDomains$n), ]

coOccuringDomains$domains <- factor(coOccuringDomains$domains, levels = coOccuringDomains$domains[order(-coOccuringDomains$n)])


#View(coOccuringDomains)


## too messy to plot all, so plot 10 most co-occuring domains

coOccuringDomains_top10 <- coOccuringDomains[c(1:10),] 


plotCoOccDoms <- ggplot(coOccuringDomains_top10, aes(x = domains, y = n, width=0.70)) +
  geom_bar(stat = "identity") + 
  labs(y = "Count", x = NULL) + 
  #scale_fill_manual(values = c("#969696","#d9d9d9")) + 
  ylim(0,1000) +
  theme_bw()
plotCoOccDoms + theme(axis.text.x = element_text(angle = 45, hjust = 1))


coOccuringDomains_top5 <- coOccuringDomains[c(1:5),] 


plotCoOccDoms_5 <- ggplot(coOccuringDomains_top5, aes(x = domains, y = n, width=0.70)) +
  geom_bar(stat = "identity") + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='identity', aes(label=n), vjust=-1) +
  #scale_fill_manual(values = c("#969696","#d9d9d9")) + 
  ylim(0,3000) +
  theme_bw()
plotCoOccDoms_5 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 8.71 x 6.42 



# look at the distributions of the top5

coOccuringDomains_top5_noRNaseH

tophit1_nRH <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR011333" & TEViralDomainID == "IPR041588")
#SKP1/BTB/POZ domain superfamily + Integrase zinc-binding domain

tophit2_nRH <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR036770" & TEViralDomainID == "IPR001995")
#Ankyrin repeat-containing domain superfamily + Peptidase A2A, retrovirus, catalytic

tophit3_nRH <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR001810" & TEViralDomainID == "IPR041426")
#F-box domain + Mos1 transposase, HTH domain

tophit4_nRH <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR011009" & TEViralDomainID == "IPR000477")
#Protein kinase-like domain superfamily + Reverse transcriptase domain

tophit5_nRH <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR036770" & TEViralDomainID == "IPR025314")
#Ankyrin repeat-containing domain superfamily + Domain of unknown function DUF4219

tophit7_nRH <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR011009" & TEViralDomainID == "IPR001584")
#Protein kinase-like domain superfamily + Integrase, catalytic core




tophit1_genus_nRH <- tophit1_nRH %>% count(genus)
tophit1_genus_nRH

tophit1_phylum_nRH <- tophit1_nRH %>% count(phylum)
tophit1_phylum_nRH
# tophit1 is mostly in vertebrates
#SKP1/BTB/POZ domain superfamily + Integrase zinc-binding domain

tophit2_genus_nRH <- tophit2_nRH %>% count(genus)
tophit2_genus_nRH

tophit2_phylum_nRH <- tophit2_nRH %>% count(phylum)
tophit2_phylum_nRH
# tophit2 is mostly in fungi
#Ankyrin repeat-containing domain superfamily + Peptidase A2A, retrovirus, catalytic

tophit3_genus_nRH <- tophit3_nRH %>% count(genus)
tophit3_genus_nRH

tophit3_phylum_nRH <- tophit3_nRH %>% count(phylum)
tophit3_phylum_nRH
# tophit3 is Caenorhabditis specific
#F-box domain + Mos1 transposase, HTH domain

tophit4_genus_nRH <- tophit4_nRH %>% count(genus)
tophit4_genus_nRH

tophit4_phylum_nRH <- tophit4_nRH %>% count(phylum)
tophit4_phylum_nRH
# tophit4 is all over the place but mostly in plants
#Protein kinase-like domain superfamily + Reverse transcriptase domain

tophit5_genus_nRH <- tophit5_nRH %>% count(genus)
tophit5_genus_nRH

tophit5_phylum_nRH <- tophit5_nRH %>% count(phylum)
tophit5_phylum_nRH
# tophit5 is only in plants, streptophyta, cool!
#Ankyrin repeat-containing domain superfamily + Domain of unknown function DUF4219


tophit7_genus_nRH <- tophit7_nRH %>% count(genus)
tophit7_genus_nRH

tophit7_phylum_nRH <- tophit7_nRH %>% count(phylum)
tophit7_phylum_nRH

## make wordclouds with top hits' genus

wordcloud2(tophit1_genus_nRH, size = 0.4)

wordcloud2(tophit2_genus_nRH, size = 0.5)

wordcloud2(tophit3_genus_nRH, size = 0.5)

wordcloud2(tophit4_genus_nRH, size = 0.5)

wordcloud2(tophit5_genus_nRH, size = 0.6)

wordcloud2(tophit7_genus_nRH, size = 0.7)

# saved as 836 x 616



## find out most promiscuous domains. as in which ubiquitous domains associate with more TE domains and vice versa.

# most promiscuous ubiquitous domains

coOccuringDomains
coOccuringDomains_mostUbiquitous <- coOccuringDomains %>% count(ubiquitousDomainID)
coOccuringDomains_mostUbiquitous


coOccuringDomains_mostUbiquitous <- coOccuringDomains_mostUbiquitous[order(-coOccuringDomains_mostUbiquitous$n), ]

coOccuringDomains_mostUbiquitous$ubiquitousDomainID <- factor(coOccuringDomains_mostUbiquitous$ubiquitousDomainID, levels = coOccuringDomains_mostUbiquitous$ubiquitousDomainID[order(-coOccuringDomains_mostUbiquitous$n)])
coOccuringDomains_mostUbiquitous

coOccuringDomains_mostUbiquitous_top10 <- coOccuringDomains_mostUbiquitous[c(1:10),]


plotCoUBIQDoms_10 <- ggplot(coOccuringDomains_mostUbiquitous_top10, aes(x = ubiquitousDomainID, y = n, width=0.70)) +
  geom_bar(stat = "identity") + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='identity', aes(label=n), vjust=-1) +
  #scale_fill_manual(values = c("#969696","#d9d9d9")) + 
  ylim(0,100) +
  theme_bw()
plotCoUBIQDoms_10 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 8.71 x 6.42 


# most promiscuous TE domains

coOccuringDomains_mostTEdoms <- coOccuringDomains %>% count(TEViralDomainID)
coOccuringDomains_mostTEdoms

coOccuringDomains_mostTEdoms <- coOccuringDomains_mostTEdoms[order(-coOccuringDomains_mostTEdoms$n), ]

coOccuringDomains_mostTEdoms$TEViralDomainID <- factor(coOccuringDomains_mostTEdoms$TEViralDomainID, levels = coOccuringDomains_mostTEdoms$TEViralDomainID[order(-coOccuringDomains_mostTEdoms$n)])
coOccuringDomains_mostTEdoms

coOccuringDomains_mostTE_top10 <- coOccuringDomains_mostTEdoms[c(1:10),]


plotCoOcTEDoms_10 <- ggplot(coOccuringDomains_mostTE_top10, aes(x = TEViralDomainID, y = n, width=0.70)) +
  geom_bar(stat = "identity") + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='identity', aes(label=n), vjust=-1) +
  #scale_fill_manual(values = c("#969696","#d9d9d9")) + 
  ylim(0,40) +
  theme_bw()
plotCoOcTEDoms_10 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 8.71 x 6.42 




##### F-box ----

AllDoms_withTaxFiltered_final_noRNaseH_singleProteins
AllDoms_withTaxFiltered_final_noRNaseH

#View(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins)
AllDoms_withTaxFiltered_final_noRNaseH_singleProteins_Fbox <- filter(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins, ubiquitousDomainID == "IPR001810")
length(unique(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins_Fbox$UniprotID))

length(unique(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins$UniprotID))
nrow(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins)

AllDoms_withTaxFiltered_final_noRNaseH_Fbox <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR001810")

length(unique(AllDoms_withTaxFiltered_final_noRNaseH$UniprotID))
nrow(AllDoms_withTaxFiltered_final_noRNaseH)

length(unique(AllDoms_withTaxFiltered_final_noRNaseH_Fbox$UniprotID))
nrow(AllDoms_withTaxFiltered_final_noRNaseH_Fbox)

plotSingleProteins_Fbox <- ggplot(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins_Fbox, aes(x = Kingdom, width=0.70, fill = Kingdom)) +
  geom_bar(aes(y = (..count..))) + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='count', aes(label=..count..), hjust=-1) +
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#54278f")) + 
  ylim(0,700) +
  coord_flip() +
  theme_bw()
plotSingleProteins_Fbox

# 8.71 x 6.42 

# how many in Caenorhabditis?
nrow(filter(AllDoms_withTaxFiltered_final_noRNaseH_singleProteins_Fbox, TEViralDomainID == "IPR041426" & genus == "Caenorhabditis"))
# 492


### to get the word clouds, chose not to do it with single proteins, because then the duplicates have different co-occurring TE domains which wouldnt be represented. 


fungi_Fbox <- filter(AllDoms_withTaxFiltered_final_noRNaseH_Fbox, Kingdom == "Fungi")

fungi2_Fbox <- fungi_Fbox %>% count(TEviralDomain)
#View(fungi2_Fbox)


top_proteins_fungi2_Fbox <- fungi2_Fbox %>%
  arrange(desc(n)) %>%  
  slice(1:10)           

wordcloud2(top_proteins_fungi2_Fbox, size = 0.2, ellipticity = 0.3)

# saved as 836 x 616

metazoa_Fbox <- filter(AllDoms_withTaxFiltered_final_noRNaseH_Fbox, Kingdom == "Metazoa")

metazoa2_Fbox <- metazoa_Fbox %>% count(TEviralDomain)
#View(metazoa2_Fbox)

top_proteins_metazoa2_Fbox <- metazoa2_Fbox %>%
  arrange(desc(n)) %>%  
  slice(1:10)         

wordcloud2(top_proteins_metazoa2_Fbox, size = 4, ellipticity = 0.3)

# saved as 836 x 616


viridiplantae_Fbox <- filter(AllDoms_withTaxFiltered_final_noRNaseH_Fbox, Kingdom == "Viridiplantae")

viridiplantae2_Fbox <- viridiplantae_Fbox %>% count(TEviralDomain)
#View(viridiplantae2_Fbox)

top_proteins_viridiplantae2_Fbox <- viridiplantae2_Fbox %>%
  arrange(desc(n)) %>%  
  slice(1:10)       

wordcloud2(top_proteins_viridiplantae2_Fbox, size = 0.12, ellipticity = 0.3)

# saved as 836 x 616

AllDoms_withTaxFiltered_final_noRNaseH_Fbox <- filter(AllDoms_withTaxFiltered_final_noRNaseH, ubiquitousDomainID == "IPR001810")

#View(AllDoms_withTaxFiltered_final_noRNaseH_Fbox)
nrow(AllDoms_withTaxFiltered_final_noRNaseH_Fbox)



plotNumberDomainsCaptured_highlight_fbox <- ggplot(AllDoms_withTaxFiltered_final_noRNaseH_Fbox, aes(x = Kingdom, width=0.70, fill = Kingdom)) +
  geom_bar(aes(y = (..count..))) + 
  labs(y = "Count", x = NULL) + 
  geom_text(stat='count', aes(label=..count..), hjust=-1) +
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#54278f")) + 
  ylim(0,700) +
  coord_flip() +
  theme_bw()
plotNumberDomainsCaptured_highlight_fbox


# 8.71 x 6.42 







