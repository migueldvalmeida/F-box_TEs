library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
string1 <- as.character(getwd())
string1 <- sub(".*/", "", string1)
pfam <- unique(read.table("IPR.tab.PFAM", sep="\t", quote=""))
colnames(pfam) <- c("protein_ID", "IPR", "description")
a <- read.table("TE_per_IPR", header = FALSE, quote="")
b <- a[,-1]
rownames(b) <- a[,1]
b[1:ncol(b)] <- lapply(b[1:ncol(b)], as.numeric)
b <- b %>% filter(rowSums(across(where(is.numeric)))!=0)
random = data.frame(t(b))
random <- random[-1,]
obs = data.frame(t(b))
obs <- obs[1,]
bobs <- data.frame(b$V2)
b$V2 <- NULL #Remove 'obs' col
c <- t(apply(b[-1], 1, function(x) sort(x, decreasing = F)))
c <- cbind(bobs, c)
d <- data.frame(t(c))
e <- data.frame(1-(apply(d[, 1:ncol(d)], 2, function(c) ecdf(c)(c))))
f <- data.frame(t(e[1,]))   
f$pfam <- rownames(f) 
f$norm <- f$X1*(length(unique(f$pfam)))
#Thresholds are based on p-values corrected for multiple comparisons
threshold1.1 = unlist(f %>% filter(norm <= (max(norm)/100*2.5)) %>% summarize(max(X1)))
threshold1.2 = unlist(f %>% filter(norm <= (max(norm)/100*97.5)) %>% summarize(max(X1)))
plot1 <- f 
plot1$species <- string1
plot1$variable <- string1

res <- plot1 %>%                    
  select(pfam, species, X1, norm) %>%
  rename(IPR = pfam,
         pvalue = X1) %>%
  inner_join(pfam, by = "IPR") %>%
  distinct()

filename <- paste0(string1,"_results", ".txt")
write.table(res, file=filename, quote=FALSE, sep='\t', col.names = F, row.names = FALSE)

p1 <- ggplot(plot1, aes(x=X1, color=species)) +
  geom_histogram(fill="white", alpha=0.01, position="identity", bins = 100) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  geom_vline(data=filter(plot1, variable=="A" | species==string1), aes(xintercept=threshold1.1), colour="orange1", linetype="dashed") + 
  geom_vline(data=filter(plot1, variable=="A"| species==string1), aes(xintercept=threshold1.2), colour="orange1", linetype="dashed") + 
  facet_wrap(~variable) +
  #ylim(0, 1000) +
  xlab("Score") + 
  ylab("Number of IPRs") +
  theme(legend.position = "none")

pdfname <- paste0(string1, "_distrib", ".pdf")
pdf(file=pdfname, width = 2.5, height = 2)
print(p1)
dev.off()

sig.PFAM <- res %>%                    
  filter(norm <= (max(norm)/100*2.5)) %>%
  select(IPR, species, pvalue, norm, description) %>%
  distinct() %>%
  arrange(pvalue)
  
filename <- paste0(string1, "_TEasso", ".txt")
write.table(sig.PFAM, file=filename, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

sig.PFAM <- res %>%                    
  filter(norm >= (max(norm)/100*97.5)) %>%
  select(IPR, species, pvalue, norm, description) %>%
  distinct() %>%
  arrange(pvalue)

filename <- paste0(string1, "_TEdep", ".txt")
write.table(sig.PFAM, file=filename, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

