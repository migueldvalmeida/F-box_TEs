library(tidyverse)
library(UpSetR)



#######################################################################
####### Load lists of TEs and make them vectors ---- ##################
#######################################################################

### first for TEs

YA_25oC_TE <- read.table("../fbxa-215_09_2021_YA/Analysis_25C/deregulatedTEs_LIST_IDS_inmutants25oC_0.05padj.txt", header = FALSE, sep = "\t")
YA_25oC_TE

YA_20oC_TE <- read.table("../fbxa-215_09_2021_YA/Analysis_20C/deregulatedTEs_LIST_IDS_inmutants20oC_0.05padj.txt", header = FALSE, sep = "\t")
YA_20oC_TE


Emb_25oC_TE <- read.table("../fbxa-215_2023_Emb/Analysis_25C/deregulatedTEs_LIST_IDS_inmutants25oC_0.05padj.txt", header = FALSE, sep = "\t")
Emb_25oC_TE

Emb_20oC_TE <- read.table("../fbxa-215_2023_Emb/Analysis_20C/deregulatedTEs_LIST_IDS_inmutants20oC_0.05padj.txt", header = FALSE, sep = "\t")
Emb_20oC_TE


YA_25oC_TE_list <- YA_25oC_TE$V1

YA_20oC_TE_list <- YA_20oC_TE$V1

Emb_25oC_TE_list <- Emb_25oC_TE$V1

Emb_20oC_TE_list <- Emb_20oC_TE$V1


listInput_TEs <- list("YA_25oC" = YA_25oC_TE_list, "YA_20oC" = YA_20oC_TE_list, "Emb_25oC" = Emb_25oC_TE_list,
                  "Emb_20oC" = Emb_20oC_TE_list)

#View(listInput_TEs)




upset(fromList(listInput_TEs), order.by = "freq")

# 6.81 x 5.31

Upset_res <- upset(fromList(listInput_TEs), order.by = "freq")
Upset_res
str(Upset_res)
Upset_res$New_data



#######################################################################
####### Load lists of genes and make them vectors ---- ################
#######################################################################


YA_25oC_up <- read.table("../fbxa-215_09_2021_YA/Analysis_25C/upregulatedgenes_inmutants25oC_0.05padj.csv", header = TRUE, sep = ",")
YA_25oC_up

YA_25oC_down <- read.table("../fbxa-215_09_2021_YA/Analysis_25C/downregulatedgenes_inmutants25oC_0.05padj.csv", header = TRUE, sep = ",")
YA_25oC_down

YA_20oC_up <- read.table("../fbxa-215_09_2021_YA/Analysis_20C/upregulatedgenes_inmutants20oC_0.05padj.csv", header = TRUE, sep = ",")
YA_20oC_up

YA_20oC_down <- read.table("../fbxa-215_09_2021_YA/Analysis_20C/downregulatedgenes_inmutants20oC_0.05padj.csv", header = TRUE, sep = ",")
YA_20oC_down


Emb_25oC_up <- read.table("../fbxa-215_2023_Emb/Analysis_25C/upregulatedgenes_inmutants25oC_0.05padj.csv", header = TRUE, sep = ",")
Emb_25oC_up

Emb_25oC_down <- read.table("../fbxa-215_2023_Emb/Analysis_25C/downregulatedgenes_inmutants25oC_0.05padj.csv", header = TRUE, sep = ",")
Emb_25oC_down

Emb_20oC_up <- read.table("../fbxa-215_2023_Emb/Analysis_20C/upregulatedgenes_inmutants20oC_0.05padj.csv", header = TRUE, sep = ",")
Emb_20oC_up

Emb_20oC_down <- read.table("../fbxa-215_2023_Emb/Analysis_20C/downregulatedgenes_inmutants20oC_0.05padj.csv", header = TRUE, sep = ",")
Emb_20oC_down


YA_25oC_up_list <- YA_25oC_up$wormgenenames
YA_25oC_down_list <- YA_25oC_down$wormgenenames

YA_20oC_up_list <- YA_20oC_up$wormgenenames
YA_20oC_down_list <- YA_20oC_down$wormgenenames

Emb_25oC_up_list <- Emb_25oC_up$wormgenenames
Emb_25oC_down_list <- Emb_25oC_down$wormgenenames

Emb_20oC_up_list <- Emb_20oC_up$wormgenenames
Emb_20oC_down_list <- Emb_20oC_down$wormgenenames


listInput_YA <- list("YA_25oC_up" = YA_25oC_up_list, "YA_25oC_down" = YA_25oC_down_list, "YA_20oC_up" = YA_20oC_up_list,
                      "YA_20oC_down" = YA_20oC_down_list)

listInput_Emb <- list("Emb_25oC_up" = Emb_25oC_up_list, "Emb_25oC_down" = Emb_25oC_down_list, "Emb_20oC_up" = Emb_20oC_up_list,
                     "Emb_20oC_down" = Emb_20oC_down_list)

listInput_ALLgenes <- list("YA_25oC_up" = YA_25oC_up_list, "YA_25oC_down" = YA_25oC_down_list, "YA_20oC_up" = YA_20oC_up_list,
                           "YA_20oC_down" = YA_20oC_down_list,"Emb_25oC_up" = Emb_25oC_up_list, "Emb_25oC_down" = Emb_25oC_down_list, "Emb_20oC_up" = Emb_20oC_up_list,
                           "Emb_20oC_down" = Emb_20oC_down_list)


upset(fromList(listInput_YA), order.by = "freq")

upset(fromList(listInput_Emb), order.by = "freq")

upset(fromList(listInput_ALLgenes), order.by = "freq", nsets = 8)




