# DAPC script
# Jamie Hudson
# Created: 12 Sep 2019
# Edited: 06 Feb 2021

library(adegenet)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)

# Set wd ------------------------------------------------------------------

setwd("~/OneDrive - University of Southampton/PhD/GBS/adegenet/dapc/pyura_praeputialis/")

# Input data --------------------------------------------------------------

pyura_neutral <- read.genepop("05Feb21_pyura_neutral.gen", ncode = 3)

pyura_neutral$pop <- as.factor(c(rep("C4",9),rep("C6",6),rep("A2",17),rep("C2",12),rep("A5",14),
                                 rep("C3",8),rep("A3",19),rep("C5",13),rep("A4",13),
                                 rep("C1",7),rep("A7",13),rep("A6",17),rep("A1",16)))

SAUS_pops <- c("A4", "A5", "A6", "A7")

pyura_neutral_noSAUS <- pyura_neutral[!pop(pyura_neutral) %in% SAUS_pops] 
pop(pyura_neutral_noSAUS)

EAUS_pops <- c("A1", "A2", "A3", "A7")

pyura_neutral_anto <- pyura_neutral_noSAUS[!pop(pyura_neutral_noSAUS) %in% EAUS_pops]
pop(pyura_neutral_anto)

# xvalDapc to calculate number of PCs to retain ---------------------------

# neutral- all sites

set.seed(999)
pyurax <- xvalDapc(tab(pyura_neutral, NA.method= "mean"), pop(pyura_neutral))
pyurax[2:6] # 40

set.seed(999)
system.time(pyurax <- xvalDapc(tab(pyura_neutral, NA.method= "mean"), pop(pyura_neutral),
                               n.pca = 35:45, n.rep = 1000,
                               parallel = "multicore", ncpus =4L))
pyurax[2:6] # 37

boxplot(pyurax$`Cross-Validation Results`$success~pyurax$`Cross-Validation Results`$n.pca, xlab = "No. PCA comp", ylab = "Class success")

# DAPC without a priori ---------------------------------------------------

#Don't know why this doesn't work
grp <- find.clusters(pyura_neutral, n.pca=300) #300 pcs and 2 clusters
dapc_pyura <- dapc(pyura_neutral, grp$grp, n.pca=37)

scatter(dapc_pyura, scree.da=FALSE, col = c("purple","orange"), legend = T)
compoplot(dapc_pyura, lab = rownames(dapc_pyura))

100*dapc_pyura$eig[1]/sum(dapc_pyura$eig)
#Population grouping

table(pop(pyura_neutral))
table(pop(pyura_neutral),grp$grp)

# DAPC with a priori ------------------------------------------------------

dapc.pyura <- dapc(pyura_neutral, n.pca=37, n.da = 10)
pdf(file = "./figures/06Feb21_dapc_pyura_neutral.pdf",
    width = 10,
    height = 6)
(A <- scatter(dapc.pyura, scree.da=FALSE, col = c("#762a83",
                                            "#762a83",
                                            "#762a83",
                                            "#1b7837",
                                            "#1b7837",
                                            "#1b7837",
                                            "#1b7837",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012"),
        cex=3, pch=20, cell=1.5, clab = 1.2))
dev.off()
compoplot(dapc.pyura, col = c("#762a83",
                              "#762a83",
                              "#762a83",
                              "#1b7837",
                              "#1b7837",
                              "#1b7837",
                              "#1b7837",
                              "#d86012",
                              "#d86012",
                              "#d86012",
                              "#d86012",
                              "#d86012",
                              "#d86012")) ## Note order of bars is same as pop(pyura_neutral)
100*dapc.pyura$eig[1]/sum(dapc.pyura$eig)

# Create plot using individual labels- note these aren't labels used in manucsript- to do this.
df <- data.frame(x = dapc.pyura$ind.coord[,1], y = dapc.pyura$ind.coord[,2])
# use the text function to add labels to the positions given by the coordinates you used in plot
s.label(dfxy = df, xax=1, yax=2,
        clabel=0.7, # change the size of the labels
        grid=T, addaxes=T) # do not draw lines or axes in addition

## Figure S2C - neutral loci, sites = Antofagasta and East Australia

set.seed(999)
pyura_noSAUSx <- xvalDapc(tab(pyura_neutral_noSAUS, NA.method= "mean"), pop(pyura_neutral_noSAUS))
pyura_noSAUSx[2:6] # 10

set.seed(999)
system.time(pyura_noSAUSx <- xvalDapc(tab(pyura_neutral_noSAUS, NA.method= "mean"), pop(pyura_neutral_noSAUS),
                               n.pca = 5:15, n.rep = 1000,
                               parallel = "multicore", ncpus =4L))
pyura_noSAUSx[2:6] # 5

dapc.pyura.noSAUS <- dapc(pyura_neutral_noSAUS, n.pca=5, n.da = 10)

pdf(file = "./figures/06Feb21_dapc_pyura_noSAUS_neutral.pdf",
    width = 10,
    height = 6)
(B <- scatter(dapc.pyura.noSAUS, scree.da=FALSE, col = c("#762a83",
                                            "#762a83",
                                            "#762a83",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012",
                                            "#d86012"),
        cex=3, pch=20, cell=1.5, clabel=1.2))
dev.off()

## Figure S2E - neutral loci, sites = Antofagasta and East Australia

set.seed(999)
pyura_noAUSx <- xvalDapc(tab(pyura_neutral_anto, NA.method= "mean"), pop(pyura_neutral_anto))
pyura_noAUSx[2:6] # 5

set.seed(999)
system.time(pyura_noAUSx <- xvalDapc(tab(pyura_neutral_anto, NA.method= "mean"), pop(pyura_neutral_anto),
                                      n.pca = 0:10, n.rep = 1000,
                                      parallel = "multicore", ncpus =4L))
pyura_noAUSx[2:6] # 6

dapc.pyura.anto <- dapc(pyura_neutral_anto, n.pca=6, n.da = 10)
pdf(file = "./figures/06Feb21_dapc_pyura_anto_neutral.pdf",
    width = 10,
    height = 6)
(C <- scatter(dapc.pyura.anto, scree.da=FALSE, col = c("#d86012",
                                                   "#d86012",
                                                   "#d86012",
                                                   "#d86012",
                                                   "#d86012"),
        cex=3, pch=20, cell=1.5, clabel=1.2))
dev.off()

# candidate_loci ----------------------------------------------------------

pyura_candidate <- read.genepop("06Feb21_pyura_candidate.gen", ncode = 3)

pyura_candidate$pop <- as.factor(c(rep("C4",9),rep("C6",6),rep("A2",17),rep("C2",12),rep("A5",14),
                                 rep("C3",8),rep("A3",19),rep("C5",13),rep("A4",13),
                                 rep("C1",7),rep("A7",13),rep("A6",17),rep("A1",16)))

set.seed(999)
pyura_candidate.x <- xvalDapc(tab(pyura_candidate, NA.method= "mean"), pop(pyura_candidate))
pyura_candidate.x[2:6] #25

set.seed(999)
system.time(pyura_candidate.x <- xvalDapc(tab(pyura_candidate, NA.method= "mean"), pop(pyura_candidate),
                               n.pca = 20:30, n.rep = 1000,
                               parallel = "multicore", ncpus =4L))
pyura_candidate.x[2:6] # 20

dapc.pyura.candidate <- dapc(pyura_candidate, n.pca=20, n.da = 10)

pdf(file = "./figures/06Feb21_dapc_pyura_candidate.pdf",
    width = 10,
    height = 6)
(D <- scatter(dapc.pyura.candidate, scree.da=FALSE, col = c("#762a83",
                                                      "#762a83",
                                                      "#762a83",
                                                      "#1b7837",
                                                      "#1b7837",
                                                      "#1b7837",
                                                      "#1b7837",
                                                      "#d86012",
                                                      "#d86012",
                                                      "#d86012",
                                                      "#d86012",
                                                      "#d86012",
                                                      "#d86012"),
        cex=3, pch=20, cell=1.5, clabel=1.2))
dev.off()
compoplot(dapc.pyura.candidate, col = c("#762a83",
                                        "#762a83",
                                        "#762a83",
                                        "#1b7837",
                                        "#1b7837",
                                        "#1b7837",
                                        "#1b7837",
                                        "#d86012",
                                        "#d86012",
                                        "#d86012",
                                        "#d86012",
                                        "#d86012",
                                        "#d86012")) ## Note order of bars is same as pop(pyura_neutral)
100*dapc.pyura$eig[1]/sum(dapc.pyura$eig)

SAUS_pops <- c("A4", "A5", "A6", "A7")

pyura_candidate_noSAUS <- pyura_candidate[!pop(pyura_candidate) %in% SAUS_pops] 
pop(pyura_candidate_noSAUS)

EAUS_pops <- c("A1", "A2", "A3", "A7")

pyura_candidate_anto <- pyura_candidate_noSAUS[!pop(pyura_candidate_noSAUS) %in% EAUS_pops]
pop(pyura_candidate_anto)

## Figure S2D - candidate loci, sites = Antofagasta and East Australia

set.seed(999)
pyura_noSAUS.candidate.x <- xvalDapc(tab(pyura_candidate_noSAUS, NA.method= "mean"), pop(pyura_candidate_noSAUS))
pyura_noSAUS.candidate.x[2:6] # 15

set.seed(999)
system.time(pyura_noSAUS.candidate.x <- xvalDapc(tab(pyura_candidate_noSAUS, NA.method= "mean"), pop(pyura_candidate_noSAUS),
                                      n.pca = 10:20, n.rep = 1000,
                                      parallel = "multicore", ncpus =4L))
pyura_noSAUS.candidate.x[2:6] # 15

dapc.pyura.candidate.noSAUS <- dapc(pyura_candidate_noSAUS, n.pca=15, n.da = 10)

pdf(file = "./figures/06Feb21_dapc_pyura_noSAUS_candidate.pdf",
    width = 10,
    height = 6)
(E <- scatter(dapc.pyura.candidate.noSAUS, scree.da=FALSE, col = c("#762a83",
                                                             "#762a83",
                                                             "#762a83",
                                                             "#d86012",
                                                             "#d86012",
                                                             "#d86012",
                                                             "#d86012",
                                                             "#d86012",
                                                             "#d86012"),
        cex=3, pch=20, cell=1.5, clabel=1.2))
dev.off()

## Figure S2F - candidate loci, sites = Antofagasta and East Australia

set.seed(999)
pyura_noAUS.candidate.x <- xvalDapc(tab(pyura_candidate_anto, NA.method= "mean"), pop(pyura_candidate_anto))
pyura_noAUS.candidate.x[2:6] # 15/25

set.seed(999)
system.time(pyura_noAUS.candidate.x <- xvalDapc(tab(pyura_candidate_anto, NA.method= "mean"), pop(pyura_candidate_anto),
                                     n.pca = 10:30, n.rep = 1000,
                                     parallel = "multicore", ncpus =4L))
pyura_noAUS.candidate.x[2:6] # 25

dapc.pyura.candidate.anto <- dapc(pyura_candidate_anto, n.pca=25, n.da = 10)

pdf(file = "./figures/06Feb21_dapc_pyura_anto_candidate.pdf",
    width = 10,
    height = 6)
(F <- scatter(dapc.pyura.candidate.anto, scree.da=FALSE, col = c("#d86012",
                                                   "#d86012",
                                                   "#d86012",
                                                   "#d86012",
                                                   "#d86012"),
        cex=3, pch=20, cell=1.5, clabel=1.2))
dev.off()


# RDA temp plot -----------------------------------------------------------

pyura_rda <- read.genepop("06Feb21_pyura_candidate_RDAsst.gen", ncode = 3)

pyura_rda$pop <- as.factor(c(rep("C4",9),rep("C6",6),rep("A2",17),rep("C2",12),rep("A5",14),
                                 rep("C3",8),rep("A3",19),rep("C5",13),rep("A4",13),
                                 rep("C1",7),rep("A7",13),rep("A6",17),rep("A1",16)))

SAUS_pops <- c("A4", "A5", "A6", "A7")

pyura_rda_noSAUS <- pyura_rda[!pop(pyura_rda) %in% SAUS_pops] 
pop(pyura_rda_noSAUS)

EAUS_pops <- c("A1", "A2", "A3", "A7")

pyura_rda_anto <- pyura_rda_noSAUS[!pop(pyura_rda_noSAUS) %in% EAUS_pops]
pop(pyura_rda_anto)

set.seed(999)
pyura_anto.candidateRDA.x <- xvalDapc(tab(pyura_rda_anto, NA.method= "mean"), pop(pyura_rda_anto))
pyura_anto.candidateRDA.x[2:6] # 10

set.seed(999)
system.time(pyura_anto.candidateRDA.x <- xvalDapc(tab(pyura_rda_anto, NA.method= "mean"), pop(pyura_rda_anto),
                                                n.pca = 10:30, n.rep = 1000,
                                                parallel = "multicore", ncpus =4L))
pyura_anto.candidateRDA.x[2:6] # 10

dapc.pyura.candidateRDA.anto <- dapc(pyura_rda_anto, n.pca=10, n.da = 10)

pdf(file = "./figures/06Feb21_dapc_pyura_anto_candidateRDA.pdf",
    width = 10,
    height = 6)
(G <- scatter(dapc.pyura.candidateRDA.anto, scree.da=FALSE, col = c("#993404","#cc4c02","#ec7014","#fe9929","#fec44f","#fc8d59"),
              cex=3, pch=20, cell=1.5, clabel=1.2))
dev.off()

sessionInfo()

