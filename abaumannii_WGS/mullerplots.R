### Plotting allele frequency plots and muller diagrams for A. baummannii populations in Scribner et al

setwd("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/allele_and_muller_plots")

library(reshape2)
library(ggplot2)

#read in 25 percent cutoff chart with added nje, manually filtered for implausible genes
abaum_snps <- read.csv(file = "/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/nje_25_filtered.csv", header=TRUE, row.names = 1, check.names = FALSE) 

#rename genes to more informative or concise names
rownames(abaum_snps) <- gsub("elongation factor G::ACX60_RS14040::", "fusA1 ", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("cytochrome ubiquinol oxidase subunit I::ACX60_RS06730::", "cyoB ", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("ubiquinol oxidase subunit II::ACX60_RS06735::", "cyoA ", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::", "ptsP ", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("stage II sporulation protein E/antiantisigma factor::", "", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("phosphocarrier protein HPr::ACX60_RS15180::", "HPr ", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("hypothetical protein/hypothetical protein::", "", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("hypothetical protein::", "", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("tRNAAsn::", "", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("50S ribosomal protein L7/L12::", "", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("NADHquinone oxidoreductase subunit B::", "", rownames(abaum_snps))
rownames(abaum_snps) <- gsub("ACX60_RS15010/ACX60_RS15015::intergenic", "ACX60_RS15010/15015", rownames(abaum_snps))

#switch so that sample names are row names 
t_abaum_snps <- as.data.frame(t(abaum_snps))
colnames(t_abaum_snps) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(t_abaum_snps))
#split rownames into sample components
sample <- colsplit(rownames(t_abaum_snps), "\\.", c("b_p", "day", "replicate", "xmic"))
t_abaum_snps <- cbind(t_abaum_snps, sample)

#subset the antibiotic samples only
t_abaum_snps <- subset(t_abaum_snps, t_abaum_snps$xmic == "inc")
#subset samples into planktonic and biofilm
snps_p <- subset(t_abaum_snps, t_abaum_snps$b_p == "P")
snps_b <- subset(t_abaum_snps, t_abaum_snps$b_p == "B")

muller <- function(snps_p, rep, b_p) {
  #subset by replicate
  snps_p <- subset(snps_p, snps_p$replicate == rep)
  #change rownames to only day
  rownames(snps_p) <- snps_p$day
  #remove unnecessary columns
  snps_p$sample <- NULL
  snps_p$day <- NULL
  snps_p$b_p <- NULL
  snps_p$replicate <- NULL
  snps_p$xmic <- NULL
  #transpose 
  snps_p <- t(snps_p)
  #convert to data.frame
  snps_p <- as.data.frame(snps_p, header=TRUE)
  #add day 0
  snps_p$"0" <- 0.0
  #convert frequencies to numeric
  snps_p[,1:ncol(snps_p)] <- apply(snps_p[,1:ncol(snps_p)], 2, function(x) as.numeric(as.character(x)))
  #add row sums and filter the mutations that don't occur in this lineage
  snps_p$sums <- rowSums(snps_p)
  snps_p <- subset(snps_p, snps_p$sums > 0 )
  #remove row sums column
  snps_p_done <- snps_p[,-ncol(snps_p)]
  #add trajectory column
  snps_p_done$Trajectory <- c(1:length(snps_p_done$"0"))
  #add gene column
  snps_p_done$Gene <- rownames(snps_p_done)
  #write into CSV file for muller
  write.csv(snps_p_done, file = paste(b_p, rep, "mullerinput.csv", sep="_")) 
  
  #prepare for allele frequency plotting, remove trajectory column
  snps_p_done$Trajectory <- NULL
  #melt
  melt_snps_done <- (melt(snps_p_done, id.vars = "Gene", value.name = "frequency"))#plot ramp
  #convert to numeric
  melt_snps_done[,2:ncol(melt_snps_done)] <- apply(melt_snps_done[,2:ncol(melt_snps_done)], 2, function(x) as.numeric(as.character(x)))
  
  snps_plot <- ggplot(melt_snps_done,aes(x=variable,y=frequency,color=Gene))+
    theme(axis.title.x = element_text(size=48),axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 1, size=32), 
          axis.text.y = element_text(colour = "black",size=32),axis.title.y = element_text(size=48), plot.title = element_text(face="bold",size = 48,hjust=0.5), 
          legend.title = element_blank(), legend.text = element_text(size = 32), legend.position="bottom",strip.text.x = element_text(size=32), 
          strip.text.y = element_text(size=32),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) +
    geom_point(size=4)+ guides(col = guide_legend(ncol = 2)) +
    ylab("Allele Frequency") +
    xlab("Day") +
    scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
    scale_y_continuous(limits = c(0,100)) +
    ggtitle(paste(b_p, rep, sep = " ")) +    
    geom_line(size=2) 
  pdf(paste(b_p, rep, ".pdf", sep = ""), width= 12, height = 10, useDingbats=F) 
  print(snps_plot) #print our plot
  dev.off()
}

muller(snps_p, 1, "Planktonic")
muller(snps_p, 2, "Planktonic")
muller(snps_p, 3, "Planktonic")

#muller plots were created with lolipop v0.5.2 using the muller input csv files generated above

#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/Planktonic_1_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/p1 
#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/Planktonic_2_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/p2
#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/Planktonic_3_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/p3

muller(snps_b, 1, "Biofilm")
muller(snps_b, 2, "Biofilm")
muller(snps_b, 3, "Biofilm")
#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/Biofilm_1_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/b1
#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/Biofilm_2_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/b2
#python /Users/mrs/muller_diagrams/mullerplot --input /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/Biofilm_3_mullerinput.csv --output /Users/mrs/Documents/PA14_Tobi_M9/tob_paper/abaum_muller/b3


