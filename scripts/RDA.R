# RDA_script.R
# File adapted from Brenna Forester's document https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# Jamie Hudson
# Created: 19 May 2020
# Edited: 12 Feb 2021


# Load packages -----------------------------------------------

library(psych)    # Used to investigate correlations among env_data.finalictors
library(vegan)    # Used to run RDA
library(adegenet)
library(dplyr)

# Load data ---------------------------------------------------------------

# Convert vcf to PLINK and then import PLINK 

plink.pyura <- read.table("../data/pyura_rda.raw", header=TRUE, row.names=1)
dim(plink.pyura)

plink.pyura <- subset(plink.pyura, select=-c(PAT, MAT, SEX, PHENOTYPE)) # Remove columns

#  Add population and region information
plink.pyura$pop <- c(rep("C4",9),rep("C6",6),rep("A2",17),rep("C2",12),rep("A5",14),
                               rep("C3",8),rep("A3",19),rep("C5",13),rep("A4",13),
                               rep("C1",7),rep("A7",13),rep("A6",17),rep("A1",16))

plink.pyura$reg <- c(rep("CHILE", 15), rep("EAUS", 17), rep("CHILE", 12),
                     rep("SAUS", 14), rep("CHILE", 8), rep("EAUS", 19),
                     rep("CHILE", 13), rep("SAUS", 13), rep("CHILE", 7),
                     rep("SAUS", 30), rep("EAUS", 16))

plink.pyura <- plink.pyura[, c(1,1207,1208,2:1206)]

# Check for missing data
sum(is.na(plink.pyura)) # 36612 NAs in the matrix (~18% missing data)

# Impute missing data -----------------------------------------------------

# RDA can't handle missing data, so must impute NA values.

# Group samples by region
plink.pyura_pop <- plink.pyura %>% 
  group_by(pop)

# Split into pops
plink.pyura_pop_group <- group_split(plink.pyura_pop)

 # Convert each split pop into a dataframe
A1_split <- as.data.frame(plink.pyura_pop_group[[1]])
A2_split <- as.data.frame(plink.pyura_pop_group[[2]])
A3_split <- as.data.frame(plink.pyura_pop_group[[3]])
A4_split <- as.data.frame(plink.pyura_pop_group[[4]])
A5_split <- as.data.frame(plink.pyura_pop_group[[5]])
A6_split <- as.data.frame(plink.pyura_pop_group[[6]])
A7_split <- as.data.frame(plink.pyura_pop_group[[7]])
C1_split <- as.data.frame(plink.pyura_pop_group[[8]])
C2_split <- as.data.frame(plink.pyura_pop_group[[9]])
C1_2_split <- rbind(C1_split, C2_split) # As there is >= 1 locus with all NAs in C1, have to combine C1 with C2 before imputing
C3_split <- as.data.frame(plink.pyura_pop_group[[10]])
C4_split <- as.data.frame(plink.pyura_pop_group[[11]])
C5_split <- as.data.frame(plink.pyura_pop_group[[12]])
C6_split <- as.data.frame(plink.pyura_pop_group[[13]])
C5_6_split <- rbind(C5_split, C6_split) # As there is >= 1 locus with all NAs in C6, have to combine C1 with C2 before imputing

# Impute NAs with the most common allele at each site
A1_rep <- apply(A1_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
A2_rep <- apply(A2_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
A3_rep <- apply(A3_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
A4_rep <- apply(A4_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
A5_rep <- apply(A5_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
A6_rep <- apply(A6_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
A7_rep <- apply(A7_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
C1_2_rep <- apply(C1_2_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
C3_rep <- apply(C3_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
C4_rep <- apply(C4_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
C5_6_rep <- apply(C5_6_split, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

plink.pyura.imp <- rbind.data.frame(A1_rep, A2_rep, A3_rep, A4_rep, A5_rep, A6_rep, A7_rep, C1_2_rep, C3_rep, C4_rep, C5_6_rep, stringsAsFactors = F)
plink.pyura.imp <- plink.pyura.imp %>% 
  group_by(pop)

plink.pyura.imp.df <- as.data.frame(sapply(plink.pyura.imp, as.numeric))

# Add individual names as rownames
rownames.pyura <- plink.pyura.imp$IID
rownames(plink.pyura.imp.df) <- rownames.pyura

# Remove extra columns now
plink.pyura.imp.df <- subset(plink.pyura.imp.df, select=-c(IID, pop, reg))

####

env_data <- read.csv("../data/pyura_environmental_data.csv")
str(env_data) # Look at the structure of the data frame
env_data$ind <- as.character(env_data$ind)

# Confirm that genotypes and environmental data are in the same order
identical(rownames(plink.pyura.imp.df), env_data[,1]) 

# Remove collinear variables
pairs.panels(env_data[,4:18], scale=T)
env_data_1 <- subset(env_data, select=-c(sst_min, cur_vel_mean, cur_vel_max, cur_vel_min))
pairs.panels(env_data_1[,4:14], scale=T)
env_data_2 <- subset(env_data_1, select=-c(chl_mean, diss_ox_min))
pairs.panels(env_data_2[,4:12], scale=T)
env_data_3 <- subset(env_data_2, select=-c(sal_max, sal_min, chl_max, diss_ox_mean))
pairs.panels(env_data_3[,4:8], scale=T)
env_data.final <- env_data_3[,4:6]
pairs.panels(env_data.final, scale = T)
str(env_data.final)

# Run RDA
pyura.rda <- rda(plink.pyura.imp.df ~ ., data=env_data.final, scale=T)
pyura.rda

RsquareAdj(pyura.rda) # Our constrained ordination explains about 2% of the variation (0.02668112); this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental env_data.finalictors (e.g., most SNPs will be neutral).

summary(eigenvals(pyura.rda, model = "constrained"))
screeplot(pyura.rda)

# Run the below code to get the significance of full model
#signif.full <- anova.cca(pyura.rda, parallel=getOption("mc.cores")) # default is permutation=999
#signif.full

# Run the below code to get the significance of axes - takes a while
#signif.axis <- anova.cca(pyura.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

# Check for Variance Inflation Factors of model - all < 5 so multicollinearity is okay
vif.cca(pyura.rda)

# Plot RDA

levels(env_data$reg) <- c("Antofagasta","East Australia","South East Australia")
reg <- env_data$reg
bg <- c("#d86012","#762a83","#1b7837") # 3 nice colors for our ecotypes
reg_1 <- gsub("anto","#d86012",reg)
reg_2 <- gsub("se_aus","#1b7837",reg_1)
reg_3 <- gsub("e_aus","#762a83",reg_2)

plot(pyura.rda, type="n", scaling=3)
points(pyura.rda, display="species", pch=4, cex=0.7, col="gray88", scaling=3)           # the SNPs
points(pyura.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=reg_3) # the squirts
text(pyura.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the env_data.finalictors
legend("bottomright", legend=levels(reg), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# Identify candidate SNPs involved in local adaptation

load.rda <- scores(pyura.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 8
cand2 <- outliers(load.rda[,2],3) # 14
cand3 <- outliers(load.rda[,3],3) # 10

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=3)  # 3 columns for 3 env_data.final predictors
colnames(foo) <- c("sst_max","sal_mean","chl_min")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- plink.pyura.imp.df[,nam]
  foo[i,] <- apply(env_data.final,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  # 2 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1

table(foo[foo[,1]==2,2]) #  no duplicates on axis 2


table(foo[foo[,1]==3,2]) # 2 duplicate on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
}

colnames(cand)[7] <- "env_data.finalictor"
colnames(cand)[8] <- "correlation"

table(cand$env_data.finalictor) 

sel <- cand$snp
env <- cand$env_data.finalictor
env[env=="sst_max"] <- '#E69F00'
env[env=="sal_mean"] <- '#56B4E9'
env[env=="chl_min"] <- '#009E73'

# color by env_data.finalictor:
col.env_data.final <- rownames(pyura.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.env_data.final)
  col.env_data.final[foo] <- env[i]
}

col.env_data.final[grep("loc",col.env_data.final)] <- '#f1eef6' # non-candidate SNPs
empty <- col.env_data.final
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#E69F00', '#56B4E9','#009E73')

# axes 1 & 2
plot(pyura.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(pyura.rda, display="species", pch=21, cex=1, col="gray32", bg=col.env_data.final, scaling=3)
points(pyura.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(pyura.rda, scaling=3, display="bp", col="#000000", cex=1)
legend("bottomright", legend=c("sst_max","sal_mean","chl_min"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(pyura.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(pyura.rda, display="species", pch=21, cex=1, col="gray32", bg=col.env_data.final, scaling=3, choices=c(1,3))
points(pyura.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(pyura.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("sst_max","sal_min","chl_min","diss_os_max","pH"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


