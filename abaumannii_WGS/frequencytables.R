#### Plotting A. baumannii mutation frequency tables for Scribner et al.

library(reshape2)
library(ggplot2)

# Read in csv file indicating mutation frequencies for SNP and New Junction Evidence variants
# File should be manually filtered for implausible genes (output of mutation_filtering.R)
abaum_snps <- read.csv(file = "/Users/mrs/Documents/tob_paper_data/abaumannii_wgs/nje_25_filtered.csv", header=TRUE, row.names = 1, check.names = FALSE) 
abaum_snps$count <- NULL
abaum_snps$sum <- NULL
abaum_snps$desc_gene_annot <- NULL
#switch so that sample names are row names 
t_abaum_snps <- as.data.frame(t(abaum_snps))
#split rownames into sample components
sample <- colsplit(rownames(t_abaum_snps), "\\.", c("b_p", "day", "replicate", "xmic"))
t_abaum_snps <- cbind(t_abaum_snps, sample)

# create day 12 frequency dataframe
day12 <- subset(t_abaum_snps, t_abaum_snps$day == "12")
#separate out antibiotic lineages
day12_ab <- subset(day12, day12$xmic == "inc")


####### Day 12, no drug and drug populations ##########

#remove sample info columns
day12$sample <- NULL
day12$day <- NULL
day12$b_p <- NULL
day12$replicate <- NULL
day12$xmic <- NULL
#convert contents of df to numeric
day12[, 1:ncol(day12)] <- sapply(day12[,1:ncol(day12)], function(x) as.numeric(as.character(x)))
#create row for sums for each mutation at day 12
day12[nrow(day12)+1, ] <- colSums(day12)
#remove genes that are at zero frequency on day 12
day12 <- day12[, which(day12[nrow(day12), ] > 0) ]
#remove colsums row
day12 <- day12[-nrow(day12), ]
#transpose
t_day12 <- as.data.frame(t(day12))
# edit gene names to be more informative
rownames(t_day12) <- gsub("ACX60_RS14040", "fusA1", rownames(t_day12))
rownames(t_day12) <- gsub(".*::fusA1", "fusA1", rownames(t_day12))
rownames(t_day12) <- gsub("ACX60_RS16050", "ptsP", rownames(t_day12))
rownames(t_day12) <- gsub(".*::ptsP", "ptsP", rownames(t_day12))
rownames(t_day12) <- gsub(".*::ACX60_RS16660", "ACX60_RS16660", rownames(t_day12))
rownames(t_day12) <- gsub(".*::ACX60_RS10565/ACX60_RS10570", "ACX60_RS10565/ACX60_RS10570", rownames(t_day12))
rownames(t_day12) <- gsub(".*::ACX60_RS06730", "cyoA", rownames(t_day12))
rownames(t_day12) <- gsub(".*::ACX60_RS06505", "ACX60_RS06505", rownames(t_day12))
ab_day12 <- as.data.frame(t_day12)


####### Day 12, drug populations only #############

#remove sample info columns
day12_ab$sample <- NULL
day12_ab$day <- NULL
day12_ab$b_p <- NULL
day12_ab$replicate <- NULL
day12_ab$xmic <- NULL
#convert contents of df to numeric
day12_ab[, 1:ncol(day12_ab)] <- sapply(day12_ab[,1:ncol(day12_ab)], function(x) as.numeric(as.character(x)))
#create row for sums for each mutation at day 12
day12_ab[nrow(day12_ab)+1, ] <- colSums(day12_ab)
#remove genes that are at zero frequency on day 12
day12_ab <- day12_ab[, which(day12_ab[nrow(day12_ab), ] > 0) ]
#remove colsums row
day12_ab <- day12_ab[-nrow(day12_ab), ]
#transpose
t_day12_ab <- t(day12_ab)
#rename columns
colnames(t_day12_ab) <- c("Bio 1", "Bio 2", "Bio 3", "Bio 4", "Bio 5", "Plank 1", "Plank 2", "Plank 3", "Plank 4", "Plank 5")
# edit gene names to be more informative, as above
rownames(t_day12_ab) <- gsub(".*::ACX60_RS14040", "fusA1", rownames(t_day12_ab))
rownames(t_day12_ab) <- gsub(".*::ACX60_RS16050", "ptsP", rownames(t_day12_ab))
rownames(t_day12_ab) <- gsub(".*::ACX60_RS16660", "ACX60_RS16660", rownames(t_day12_ab))
rownames(t_day12_ab) <- gsub(".*::ACX60_RS06730", "cyoA", rownames(t_day12_ab))
#invert order 
day12_ab <- t_day12_ab[order(nrow(t_day12_ab):1), ]

#create df for aggregating by gene
t_day12_ab_gene <- t_day12_ab
#remove mutation information from rownames so that it is only gene
rownames(t_day12_ab_gene) <- gsub("::.*", "", rownames(t_day12_ab))
#aggregate by gene
day12_ab_gene <- aggregate(t_day12_ab_gene, by=list(as.character(rownames(t_day12_ab_gene))), FUN=sum)
#make the gene name column the row names, then remove the gene name column
rownames(day12_ab_gene) <- day12_ab_gene$Group.1
day12_ab_gene$Group.1 <- NULL
# invert rows
day12_ab_gene <- day12_ab_gene[order(nrow(day12_ab_gene):1), ]

###############################################

# Combine A. baumannii and P. aeruginosa day 12 ab gene aggregated dfs to create figure2

# Read in P. aeruginosa
pa_ab_gene <- read.csv(file="/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/day12_ab_gene.csv", header=TRUE, row.names = 1, check.names = FALSE)
colnames(pa_ab_gene) <- c("plank1pa", "plank2pa", "plank3pa", "plank4pa", "plank5pa", "bio1pa", "bio2pa", "bio3pa", "bio4pa", "bio5pa")
day12_ab_gene$gene <- rownames(day12_ab_gene)
pa_ab_gene$gene <- rownames(pa_ab_gene)
all_ab_gene <- merge(day12_ab_gene, pa_ab_gene, by = "gene", all = TRUE)
all_ab_gene[is.na(all_ab_gene)] <- 0
rownames(all_ab_gene) <- all_ab_gene$gene
all_ab_gene$gene <- NULL
all_ab_gene[, 1:ncol(all_ab_gene)] <- sapply(all_ab_gene[,1:ncol(all_ab_gene)], function(x) as.numeric(as.character(x)))
all_ab_gene <- all_ab_gene[, c(6,7,8,9,10,1,2,3,4,5,11,12,13,14,15,16,17,18,19,20)]
all_ab_gene <- all_ab_gene[c(3,9,2,1,10,8,7,6,5,4), ]
all_ab_gene <- all_ab_gene[nrow(all_ab_gene):1, ]

all_ab_gene <- as.matrix(all_ab_gene)
all_ab_gene <- melt(all_ab_gene, id.vars = c("population", "gene"), value.name="frequency")
day12_ggplot <- ggplot(all_ab_gene, aes(Var2, Var1)) + geom_tile(aes(fill = frequency), colour = "white") + 
  scale_fill_gradient(low = "white", high = "black", limits=c(0,110)) + theme_classic() +
  theme(axis.text.x=element_text(angle = 90), axis.text=element_text(size=10, face="bold.italic"), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        panel.border = element_rect(color = "black", size = 1, fill = NA))
print(day12_ggplot)
ggsave("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/Figure2.pdf", plot = (day12_ggplot), height = 4, width = 8, device = "pdf")

###########################

#Combine A. baumannii and P. aeruginosa day 12 dfs to create supplemental figure 2
pa_day12 <- read.csv(file="/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/day12.csv", header=TRUE, row.names = 1, check.names = FALSE)
colnames(pa_day12) <- c("plank1pa notob", "plank2pa notob", "plank3pa notob","bio1pa notob", "bio2pa notob", "bio3pa notob", "plank1pa", "plank2pa", "plank3pa", "plank4pa", "plank5pa", "bio1pa", "bio2pa", "bio3pa", "bio4pa", "bio5pa")
pa_day12$gene <- rownames(pa_day12)

#change rownames to just gene_annot
gene_annot <- colsplit(rownames(pa_day12), "::", c("gene", "annot", "desc"))
gene_annot$gene_annot <- paste(gene_annot$gene, gene_annot$annot, sep=" ")
pa_day12$gene <- gene_annot$gene_annot

#add a gene column to the abaum data from above
ab_day12$gene <- rownames(ab_day12)
gene_annot <- colsplit(rownames(ab_day12), "::", c("gene", "annot"))
gene_annot$gene_annot <- paste(gene_annot$gene, gene_annot$annot, sep=" ")
ab_day12$gene <- gene_annot$gene_annot

all_day12 <- merge(ab_day12, pa_day12, by = "gene", all = TRUE)
all_day12[is.na(all_day12)] <- 0
rownames(all_day12) <- all_day12$gene
all_day12$gene <- NULL
all_day12[, 1:ncol(all_day12)] <- sapply(all_day12[,1:ncol(all_day12)], function(x) as.numeric(as.character(x)))
all_day12 <- all_day12[, c(14,15,16,11,12,13,6,7,8,9,10,1,2,3,4,5,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)]
grep("orf", rownames(all_day12))
all_day12 <- all_day12[c(5:20,34:57,4,3,29,23:27,1,2,21,22,28,30:33,59),]
all_day12 <- all_day12[nrow(all_day12):1, ]

all_day12 <- as.matrix(all_day12)
write.csv(all_day12, "/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/FigureS2.csv")
all_day12 <- melt(all_day12, id.vars = c("population", "gene"), value.name="frequency")
day12_ggplot <- ggplot(all_day12, aes(Var2, Var1)) + geom_tile(aes(fill = frequency), colour = "white") + 
  scale_fill_gradient(low = "white", high = "black", limits=c(0,100)) + theme_classic() +
  theme(axis.text.x=element_text(angle = 90), axis.text=element_text(size=12, face="bold.italic"), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        panel.border = element_rect(color = "black", size = 1, fill = NA))
print(day12_ggplot)
ggsave("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/FigureS2.pdf", plot = (day12_ggplot), height = 16, width = 16, device = "pdf")
