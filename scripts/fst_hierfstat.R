#hierfstat script to calculate pairwise Fst values
#Jamie Hudson
# Created 24 Sep 2019
# Modified 11 Feb 2021

library(adegenet)
library(hierfstat)
library(reshape2)
library(tidyverse)

# Set wd ------------------------------------------------------------------

setwd("/GBS/hierfstat/pyura_praeputialis/")


# Input data --------------------------------------------------------------

neutral.gen <- read.genepop("/05Feb21_pyura_neutral.gen", ncode = 3) # read in genepop file

neutral.hierfstat <- genind2hierfstat(neutral.gen) #convert genind to hierfstat df

# Set NAs to 0?
# neutral.hierfstat[is.na(neutral.hierfstat)] <- 0

#Read in pop names to match what is in ms
neutral.hierfstat$pop <- as.factor(c(rep("C4",9),rep("C6",6),rep("A2",17),rep("C2",12),rep("A5",14),
                                rep("C3",8),rep("A3",19),rep("C5",13),rep("A4",13),
                                rep("C1",7),rep("A7",13),rep("A6",17),rep("A1",16)))

neutral_pwfst <- genet.dist(neutral.hierfstat, method = "WC84")
str(neutral_pwfst)

neutral_pwfst_matrix<- as.matrix(neutral_pwfst)

write.table(neutral_pwfst_matrix, file="neutral_pwfst_matrix",row.names = T,col.names = T) #rename to something appropriate

# in case you are reading in from file
#neutral_pwfst_matrix <- as.matrix(read.table("allsites_neutral_pwfst_matrix"))

neutral_pwfst_matrix <- round(neutral_pwfst_matrix,4)
neutral_pwfst.melt <- melt(neutral_pwfst_matrix, na.rm =TRUE)
summary(neutral_pwfst.melt$value)

# Benjamini-Yekutieli FDR correction - a' = a/sum(1/i) where i = number of comparisons
# count number of comparisons for correction
sum(1:12) # 78 - note not sum(1:13), as despite 13 pops we don't do a pw comparison of each pop on itself
x <- seq(1,78,1) # 1 to n of number of comparisons
recip.x <- c(1/x)
sum.recip.x <- sum(recip.x)

# Bootstrap for CIs
pwfst_boot <- boot.ppfst(dat=neutral.hierfstat, quant = c((0.05/sum.recip.x),(1-(0.05/sum.recip.x))), nboot = 10000)
pwfst_boot_lci <- as.matrix(pwfst_boot[[2]])
pwfst_boot_uci <- as.matrix(pwfst_boot[[3]])

melted_pwfst_lci <- melt(pwfst_boot_lci, na.rm =F)
melted_pwfst_uci <- melt(pwfst_boot_uci, na.rm =F)

melted_pwfst_lci <- as.tibble(melted_pwfst_lci)
melted_pwfst_lci <- rename(melted_pwfst_lci, lci = value)

melted_pwfst_uci <- as.tibble(melted_pwfst_uci)
melted_pwfst_uci <- rename(melted_pwfst_uci, uci = value)

fst_ci <- cbind(melted_pwfst_uci, melted_pwfst_lci)
neutral_pwfst.melt.ci <- cbind(fst_ci, neutral_pwfst.melt)
neutral_pwfst.melt.ci <- neutral_pwfst.melt.ci[-c(1:2,4:5)]
neutral_pwfst.melt.ci <- neutral_pwfst.melt.ci[!(is.na(neutral_pwfst.melt.ci$uci) & neutral_pwfst.melt.ci$Var2 != "C6"), ] #Remove for the plot
neutral_pwfst.melt.ci$value[neutral_pwfst.melt.ci$value == 0] <- NA
nrow(neutral_pwfst.melt.ci[which(neutral_pwfst.melt.ci$uci >0 & neutral_pwfst.melt.ci$lci >0),]) # number significant comparisons

# Produce plot
(pwfst.plot <- ggplot(data = neutral_pwfst.melt.ci, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ 
    scale_fill_viridis_c(option = "magma",direction = -1, name = expression(F[ST]), na.value = "white")  +
    scale_y_discrete(limits = rev(levels(neutral_pwfst.melt$Var2))) +
    ggtitle(expression(atop("Pairwise FST, WC (1984), Neutral loci", atop(italic("1,140 loci"), ""))))+
    labs( x = "Sampling Site", y = "Sampling Site") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),axis.text.y = element_text(size = 13)) + 
    theme(axis.title = element_text(size = 16),legend.text = element_text(size =20), legend.title = element_text(size =30)) +
    theme(plot.title = element_text(size = 17)) +
    coord_fixed() +
    geom_text(aes(label = ifelse(uci >0 & lci >0, "*","")), size = 10, colour = "white", vjust = 0.75)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14),legend.text = element_text(size =10), legend.title = element_text(size =16)))
