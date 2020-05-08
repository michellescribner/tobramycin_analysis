# Code to generate frequency tables and muller plots for P. aerguionsa populations in Scribner et al

setwd("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS")

library(ggplot2)
library(ggmuller)
library(reshape2)

# Use this script to generate mutation frequency tables and Muller plots from breseq data for populations over time #

# A note on filtering: to be included in these plots, mutations must have reached a total frequency of 25% across all timepoints and lineages in the experiment. 
# Mutations in no-antibiotic lineages were including in this filtering. Mutations were filtered by position, not by gene. 
# Mutations manually determined to be biologically implausible were also removed. 
# High quality mutations reported in New Junction Evidence were included if they were related to known tobramycin resistance mechanisms.
# For P. aerugionsa, one mutation in 23S rRNA was at high frequencies and occurred in consistent populations so was determined to meet this criteria.

# Read in previously generated csv file of snps after filtering (mutation_filtering.R), with nje manually added
snps <- read.csv("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/t_snps_25_nje.csv", header=T, check.names = FALSE, row.names = 1)

#split sample column into components
sample <- colsplit(rownames(snps), ", ", c("day", "pop", "b_p", "xmic"))
sample$day <- gsub( "Day ", "", as.character(sample$day))
sample$pop <- gsub( "Pop ", "", as.character(sample$pop))
snps <- cbind(sample, snps)
#remove sums row
snps <- snps[-nrow(snps), ]

# create day 12 frequency dataframe
day12 <- subset(snps, snps$day == "12")
# subset abx samples from day 12 df
day12_ab <- subset(day12, day12$xmic == "Inc")
#subset abx samples (all days)
snps_ab <- subset(snps, snps$xmic == "Inc")

######### create dataframe showing frequency of each mutation at day 12 ###########

#remove sample info columns
day12 <- day12[,-c(1:4)]
#create row for sums for each mutation at day 12
day12[nrow(day12)+1, ] <- colSums(day12)
#remove genes that are at zero frequency on day 12
day12 <- day12[, which(day12[nrow(day12), ] > 0) ]
#remove colsums row
day12 <- day12[-nrow(day12), ]
#transpose
t_day12 <- t(day12)
#order columns
day12 <- t_day12[ ,c(4,8,12,2,6,10,3,7,11,14,16,1,5,9,13,15)]
write.csv(day12, file= "/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/day12.csv")

########### create a dataframe with only antibiotic samples at day 12 ##########

#remove sample info columns
day12_ab <- day12_ab[,-c(1:4)]
#create row for sums for each mutation at day 12
day12_ab[nrow(day12_ab)+1, ] <- colSums(day12_ab)
#remove genes that are at zero frequency on day 12
day12_ab <- day12_ab[, which(day12_ab[nrow(day12_ab), ] > 0) ]
#remove colsums row
day12_ab <- day12_ab[-nrow(day12_ab), ]
#transpose
t_day12_ab <- t(day12_ab)
#order columns
day12_ab <- t_day12_ab[ ,c(2,4,6,8,10,1,3,5,7,9)]

### AGGREGATE BY GENE
#remove mutation information from rownames so that it is only gene
day12_ab_gene <- day12_ab
rownames(day12_ab_gene) <- gsub(":.*", "", rownames(day12_ab_gene))
#aggregate by gene
day12_ab_gene <- aggregate(day12_ab_gene, by=list(rownames(day12_ab_gene)), FUN=sum)
#make the gene name column the row names, then remove the gene name column
rownames(day12_ab_gene) <- day12_ab_gene$Group.1
day12_ab_gene$Group.1 <- NULL

#put the rows in the desired order
day12_ab_gene <- day12_ab_gene[c(1,7,8,6,5,4,3,2), ]
#ggplot inverts the rows, so invert the row ahead of time so that it's right
day12_ab_gene <- day12_ab_gene[order(nrow(day12_ab_gene):1), ]
write.csv(day12_ab_gene, file= "/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/day12_ab_gene.csv")

################################

## MULLER PLOTS ##

#version conda install pygraphviz was needed
# You are using pip version 19.0.3
# Version 0.5.2 of muller scripts, updated on 8.7.19

setwd("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/allele_and_muller_plots")

snps_ab <- as.data.frame(t(snps_ab))
#change rownames to just gene_annot
gene_annot <- colsplit(rownames(snps_ab), "::", c("gene", "annot", "desc"))
gene_annot$gene_annot <- paste(gene_annot$gene, gene_annot$annot, sep=" ")
rownames(snps_ab) <- gene_annot$gene_annot
snps_ab <- as.data.frame(t(snps_ab))
colnames(snps_ab) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(snps_ab))

#subset samples into planktonic and biofilm
snps_p <- subset(snps_ab, snps_ab$b_p == "Plank")
snps_b <- subset(snps_ab, snps_ab$b_p == "Bio")

muller <- function(snps, rep, b_p) {
  #subset by replicate
  snps <- subset(snps, snps$pop == rep)
  #transpose 
  snps <- t(snps)
  #convert to data.frame
  snps <- as.data.frame(snps, header=TRUE)
  #make column names the days
  colnames(snps) <- as.character(unlist(snps[1,]))
  #add day 0
  snps$"0" <- 0.0
  #remove sample info columns
  snps_pretty <- snps[-c(1:4),]
  #convert frequencies to numeric
  snps_pretty[,1:ncol(snps_pretty)] <-(apply(snps_pretty[,1:ncol(snps_pretty)], 2, function(x) as.numeric(as.character(x)))) 
  #add row sums and filter the mutations that don't occur in this lineage
  snps_pretty$sums <- rowSums(snps_pretty)
  snps_pretty <- subset(snps_pretty, snps_pretty$sums > 0 )
  #remove row sums column
  snps_done <- snps_pretty[,-ncol(snps_pretty)]
  #add trajectory column
  snps_done$Trajectory <- (c(1:length(snps_done$"0")))
  #add gene name column
  snps_done$Gene <- rownames(snps_done)
  #write into CSV file for muller
  write.csv(snps_done, file = paste(b_p, rep, "mullerinput.csv", sep="_"))
  
  #prepare for allele frequency plotting, remove trajectory column
  snps_done$Trajectory <- NULL
  #melt
  melt_snps_done <- (melt(snps_done, id.vars = "Gene", value.name = "frequency"))#plot ramp
  #convert to numeric
  melt_snps_done[,2:ncol(melt_snps_done)] <- apply(melt_snps_done[,2:ncol(melt_snps_done)], 2, function(x) as.numeric(as.character(x)))

  snps_plot <- ggplot(melt_snps_done,aes(x=variable,y=frequency,color=Gene))+
    theme(axis.title.x = element_text(size=48),axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 1, size=32), 
    axis.text.y = element_text(colour = "black",size=32),axis.title.y = element_text(size=48), plot.title = element_text(face="bold",size = 48, hjust = 0.5), 
    legend.title = element_blank(), legend.text = element_text(size = 32), legend.position="bottom",strip.text.x = element_text(size=32), 
    strip.text.y = element_text(size=32),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) +
    geom_point(size=4)+ guides(col = guide_legend(ncol = 2)) +
    ylab("Allele Frequency") +
    xlab("Day") +
    scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
    scale_y_continuous(limits = c(0,100)) +
    ggtitle(paste(b_p, rep, sep = " ")) +    
    geom_line(size=3) 
  pdf(paste(b_p, rep, ".pdf", sep = ""), width= 12, height = 10, useDingbats=F) 
  print(snps_plot) #print our plot
  dev.off()
}

muller(snps_p, 1, "Planktonic")
muller(snps_p, 2, "Planktonic")
muller(snps_p, 3, "Planktonic")
muller(snps_b, 1, "Biofilm")
muller(snps_b, 2, "Biofilm")
muller(snps_b, 3, "Biofilm")

# Running Muller scripts on the command line:

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/Biofilm_1_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/b1

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/Biofilm_2_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/b2
# In this population, the frequencies of orfN mutations reported in Breseq's RA evidence column were slightly different from the frequencies observed by manually viewing read pileups and by clicking "RA" in breseq html format.  
# On days 10 and 12, the orfN mutation was reported in the RA column to be 100%, but clicking "RA" in breseq html format revealed that the mutation was actually at ~80-90% on these days.
# The frequencies in the RA evidence column were biologically implausible, in that these mutations were reported to be at 100%, but other mutations known by clonal WGS to be on different genetic backgrounds were also present in the population at those timepoints. 
# On the other hand, the frequencies manually calculated by examining read pileups were consistent with trajectories of other genotypes within the population, suggesting that these are more accurate. 
# We therefore adjusted the mutation frequencies from 100% to ~90% for the muller plots as these likely better reflect the allele frequencies at these timepoints. The original frequencies may still be viewed in the allele frequency plots for transparency.
# The orfN mutation for which these adjustments were required was a new junction in a repeat region, so while SNP calls appear to be reported at generally reliable frequencies for muller and allele frequency plots, junction calls were manually analyzed for accuracy of breseq's frequency prediction for use in muller plots.

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/Biofilm_3_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/b3

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/Planktonic_1_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/p1

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/Planktonic_2_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/p2
# As in biofilm population 2, the frequencies of a new junction mutation in ptsP reported in Breseq's RA evidence columnn were slightly different from frequencies observed by manually viewing read pileups and by clicking "RA" in breseq html format. 
# The mutation was called at 100% in the RA column on days 10 and 12, but clicking on the mutation revealed that it was actually ~80-90%. We adjusted these frequencies in these muller plots for the reasons indicated above. 

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/Planktonic_3_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/pa_muller/p3 
#ran by chris deitrick

