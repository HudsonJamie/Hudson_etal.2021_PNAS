# diveRsity script
# Jamie Hudson
# Created: 12 Sep 2019
# Edited: 06 Feb 2021

library(diveRsity)
library(adegenet)
library(tidyverse)

# Fis ---------------------------------------------------------------------


popgen_stats <- divBasic("../data/pyura_neutral.gen", outfile = NULL, gp = 3, 
                     bootstraps = 10000)

save(popgen_stats, file = "popgen_stats.RData") # something suitable
load("divBasic_pyura_neutral_diveRsity.RData")

###explore results stored in variable popgen_stats
names(popgen_stats)

# popgen_stats contains number of objects but interest is in Fis now.
# Fis contains a dataframe for each population with Fis value and 95% CI
# (standard and bias corrected)

#### First extract the relevant information from our results. Note BC_ refers to bias corrected confidence intervals (http://rstudio-pubs-static.s3.amazonaws.com/12475_c81fc78c66324a40b1509ee7ce3a15f6.html)
ciTable <- lapply(popgen_stats$fis, function(x){
  return(c(x["overall", "fis"], x["overall", "BC_lower_CI"],
           x["overall", "BC_upper_CI"]))
})
### convert this into a dataframe
ciTable <- as.data.frame(do.call("rbind", ciTable))
pop_names <- c("C4", "C6", "A2", "C2", "A5", "C3",
               "A3", "C5", "A4", "C1", "A7", "A6",
               "A1")
dimnames(ciTable) <- list(pop_names,
                          c("Fis", "lower", "upper"))

# inspect the table
ciTable
ciTable[order(row.names(ciTable)),]

# plot results
ggplot(ciTable, aes(x = pop_names, y = Fis)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  theme(axis.text.x=element_text(angle = -90))

# Extract Ho and He too
Ho <- as.data.frame(popgen_stats$Ho)
He <- as.data.frame(popgen_stats$He)

mean_Ho <- Ho %>% summarise_all(mean, na.rm = T)
colnames(mean_Ho) <- c("C4", "C6", "A2", "C2", "A5", "C3", 
                       "A3", "C5", "A4", "C1", "A7", "A6", 
                       "A1") # same order as input file
mean_Ho
mean_Ho[,order(colnames(mean_Ho))]

mean_He <- He %>% summarise_all(mean, na.rm = T)
colnames(mean_He) <- c("C4", "C6", "A2", "C2", "A5", "C3",
                       "A3", "C5", "A4", "C1", "A7", "A6",
                       "A1") # same order as input file
mean_He[,order(colnames(mean_He))]

# Count number of private alleles

priv_allele.pyura <- readGenepop(infile = "../adegenet/dapc/pyura_praeputialis/05Feb21_pyura_neutral.gen", gp = 3)
priv_allele.pyura$pop_names <- c("C4", "C6", "A2", "C2", "A5", "C3", 
                                 "A3", "C5", "A4", "C1", "A7", "A6", 
                                 "A1")

priv_allele.pyura$pop <- as.factor(c(rep("C4",9),rep("C6",6),rep("A2",17),rep("C2",12),rep("A5",14),
                             rep("C3",8),rep("A3",19),rep("C5",13),rep("A4",13),
                             rep("C1",7),rep("A7",13),rep("A6",17),rep("A1",16)))


allele.freq.pyura <- priv_allele.pyura$allele_freq

length(allele.freq.pyura) # number of alleles

pyura.allele.zeros <- lapply(allele.freq.pyura, function(x) which(x==0)) # Tells us which population and allele combination have a frequency of 0

pyura.private.alleles <- lapply(pyura.allele.zeros, function(x) which(length(x) == 12)) # 12 zeros means only one pop has that allele

length(which(pyura.private.alleles==1)) # Number of private alleles

pyura.private.alleles.list <- allele.freq.pyura[which(pyura.private.alleles==1)]

private.alleles.df <- data.frame(matrix(unlist(pyura.private.alleles.list), nrow=length(pyura.private.alleles.list), byrow=TRUE))

colnames(private.alleles.df) <- as.factor(c(rep("C4",2),rep("C6",2),rep("A2",2),rep("C2",2),rep("A5",2),
                                rep("C3",2),rep("A3",2),rep("C5",2),rep("A4",2),
                                rep("C1",2),rep("A7",2),rep("A6",2),rep("A1",2)))

private.alleles.df_2 <- private.alleles.df[seq(1, ncol(private.alleles.df),2)]

# For each population count number of private alleles- got to make this nicer
private.alleles.df_2 %>% 
  summarise(count = sum(A1 != 0 & A1 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(A2 != 0 & A2 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(A3 != 0 & A3 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(A4 != 0 & A4 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(A5 != 0 & A5 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(A6 != 0 & A6 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(A7 != 0 & A7 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(C1 != 0 & C1 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(C2 != 0 & C2 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(C3 != 0 & C3 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(C4 != 0 & C4 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(C5 != 0 & C5 != 1))
private.alleles.df_2 %>% 
  summarise(count = sum(C6 != 0 & C6 != 1))
                                
sessionInfo()
