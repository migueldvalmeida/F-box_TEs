library(tidyverse)
library(reshape2)
library(ggrepel)



###############################################################################
########## ---- Plotting hits as a scatterplot - after filtering  #############
###############################################################################


AllDoms_withTaxFiltered_final

Filtered_hits <- AllDoms_withTaxFiltered_final %>% count(ubiquitousDomainID)
Filtered_hits


AllUbiquitiousDoms <- read.table(file = "../OriginalDomainsToOverlap/AllDomainsKeptForTheAnalysis/AllUbiquitousDomains.txt", sep = "\t",
                      header = FALSE, fill =TRUE, quote = NULL)

#View(AllUbiquitiousDoms)
nrow(AllUbiquitiousDoms)

# remove last column
AllUbiquitiousDoms <- AllUbiquitiousDoms[,1:8]

length(unique(AllUbiquitiousDoms$V1))


nrow(AllUbiquitiousDoms)


########## Add additional taxonomic information to the dataframe 



# prepared database of taxonomizr, only need to do it once: prepareDatabase('accessionTaxa.sql')


# Extract taxonomy information from classification results
Taxonomy_AllUbiquitousDomsV4 <-getTaxonomy(AllUbiquitiousDoms$V4,'accessionTaxa.sql')

#View(Taxonomy_AllUbiquitousDomsV4)

nrow(Taxonomy_AllUbiquitousDomsV4)
nrow(AllUbiquitiousDoms)

AllUbiquitousDoms_withTax <- cbind(AllUbiquitiousDoms,Taxonomy_AllUbiquitousDomsV4)
#View(AllUbiquitousDoms_withTax)
nrow(AllUbiquitousDoms_withTax)

is.data.frame(AllUbiquitousDoms_withTax)



########## Remove metagenome and prokaryotes 


# remove prokaryotes and viruses 
unique(AllUbiquitousDoms_withTax$superkingdom)

nrow(AllUbiquitousDoms_withTax)
AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTax[!grepl("Bacteria",AllUbiquitousDoms_withTax$superkingdom),]
nrow(AllUbiquitousDoms_withTaxFiltered)
unique(AllUbiquitousDoms_withTaxFiltered$superkingdom)

AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTaxFiltered[!grepl("Archaea",AllUbiquitousDoms_withTaxFiltered$superkingdom),]
nrow(AllUbiquitousDoms_withTaxFiltered)
unique(AllUbiquitousDoms_withTaxFiltered$superkingdom)

AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTaxFiltered[!grepl("Viruses",AllUbiquitousDoms_withTaxFiltered$superkingdom),]
nrow(AllUbiquitousDoms_withTaxFiltered)
unique(AllUbiquitousDoms_withTaxFiltered$superkingdom)

#View(AllUbiquitousDoms_withTaxFiltered)

# What are the superkingdom NAs?
AllUbiDoms_withTaxFiltered_NAs <- AllUbiquitousDoms_withTaxFiltered[is.na(AllUbiquitousDoms_withTaxFiltered$superkingdom),]
#View(AllUbiDoms_withTaxFiltered_NAs)

# most of these are metagenomes... remove
nrow(AllUbiDoms_withTaxFiltered_NAs) # 55319 lines with metagenomes and uncultured organisms

AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTaxFiltered[!grepl("metagenome",AllUbiquitousDoms_withTaxFiltered$V5),]
nrow(AllUbiquitousDoms_withTaxFiltered)
unique(AllUbiquitousDoms_withTaxFiltered$superkingdom)

# removed 48759 metagenomes

AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTaxFiltered[!grepl("uncultured",AllUbiquitousDoms_withTaxFiltered$V5),]
nrow(AllUbiquitousDoms_withTaxFiltered)
unique(AllUbiquitousDoms_withTaxFiltered$superkingdom)

# removed 1441 uncultured organisms. 

AllDoms_withTaxFiltered_NAs2 <- AllUbiquitousDoms_withTaxFiltered[is.na(AllUbiquitousDoms_withTaxFiltered$superkingdom),]
#View(AllDoms_withTaxFiltered_NAs2)

# the rest are taxonomy entries with NAs in all fields, remove them all. 
unique(AllDoms_withTaxFiltered_NAs2$V5)
nrow(AllDoms_withTaxFiltered_NAs2)

# filter out all the NAs, better for simplicity, so I dont have to input data manually.

AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTaxFiltered[!is.na(AllUbiquitousDoms_withTaxFiltered$superkingdom),]
nrow(AllUbiquitousDoms_withTaxFiltered)
unique(AllUbiquitousDoms_withTaxFiltered$superkingdom)

# noticed that 1207 proteins somehow lost the protein domain entry in column V7, remove them from the dataframe.
empty_spaces <- AllUbiquitousDoms_withTaxFiltered$V7 == ""
df_with_empty_spaces <- AllUbiquitousDoms_withTaxFiltered[empty_spaces, ]
#View(df_with_empty_spaces)

AllUbiquitousDoms_withTaxFiltered <- AllUbiquitousDoms_withTaxFiltered[!(AllUbiquitousDoms_withTaxFiltered$V7 == ""), ]

nrow(AllUbiquitousDoms_withTaxFiltered)

## 11,446,044 multicellular eukaryotes proteins screened! ----

nrow(AllUbiquitousDoms_withTaxFiltered)
length(unique(AllUbiquitousDoms_withTaxFiltered$V1))

## 9,942,825 unique proteins, and 1,503,219 that have more than 1 domain in the list ----


Filtered_hits


TotalProteinsDomain <- AllUbiquitousDoms_withTaxFiltered %>% count(V7)
TotalProteinsDomain

colnames(TotalProteinsDomain) <- c("ubiquitousDomainID","n_total")
TotalProteinsDomain

forScatterPlot <- left_join(Filtered_hits,TotalProteinsDomain, by = "ubiquitousDomainID")


# now left join another dataframe to attach the class of the domain

InterproClass <- ListHits[,c(1,3)]
InterproClass

forScatterPlot <- left_join(forScatterPlot,InterproClass,by = c("ubiquitousDomainID" = "Interpro"))

# statistical model of regression

model <- lm(n ~ n_total, data = forScatterPlot)
model

summary(model)

## Adjusted R-squared:  0.7239 , p-value: 2.176e-10


ggplot(forScatterPlot, aes(x=n_total, y=n, label = Class)) +
  geom_point(size=2) + 
  geom_text_repel() + 
  labs(x = "Proteins in Interpro Domain", y = "Proteins co-occuring with TEs") +
  stat_smooth(method = lm) +
  theme_bw()




