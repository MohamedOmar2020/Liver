###############################################################################################
## Mohamed Omar
## 21/07/2019
## Goal: Discovery and validation of a small gene signature that can predict the metastasis potential in primary prostate cancer
################################################################################################

## Clean work space
rm(list = ls())

## Set the working directory
#setwd("/Volumes/Macintosh/Research/Projects/Liver")


## Load necessary packages
library(MetaIntegrator)
library(GEOquery)
library(pROC)
library(caret)
library(genefilter)
library(mltools)
library(hgu133plus2.db)
library(annotate)
library(precrec)
library(patchwork)
library(edgeR)
library(limma)

########################################################


#FattyLiverData <- getGEOData(c("GSE151158", "GSE83452", "GSE126848", "GSE130970"))
# CellLinesData <- getGEOData(c("GSE39919", "GSE70394", "GSE74577", "GSE74492"))
# 
# # Get each dataset separately
# Dataset1 <- FattyLiverData$originalData$GSE151158
# Dataset2 <- FattyLiverData$originalData$GSE83452
# Dataset3 <- FattyLiverData$originalData$GSE126848
# Dataset4 <- FattyLiverData$originalData$GSE130970
# 
# 
# ## Get the expression of dataset3
# expr3 <- getGEOSuppFiles("GSE126848")
# expr3 <- read.delim("./Data/GSE126848_Gene_counts_raw.txt")
# Dataset3$expr <- expr3
# Dataset3$keys <- expr3$key
# 
# 
# ## Get the expression of dataset4
# expr4 <- getGEOSuppFiles("GSE130970", filter_regex = "GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz")
# expr4 <- read.csv("./Data/GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv")
# Dataset4$expr <- expr4
# Dataset4$keys <- expr4$entrez_id

#save(Dataset1, Dataset2, Dataset3, Dataset4, file = "./Data/FattyLiverData.rda")



# ValDataset1 <- CellLinesData$originalData$GSE39919
# ValDataset2 <- CellLinesData$originalData$GSE70394
# ValDataset3 <- CellLinesData$originalData$GSE74577
# ValDataset4 <- CellLinesData$originalData$GSE74492
# 

#CanDataset1 <- getGEOData(c("GSE65801"))
#CanDataset1 <- CanDataset1$originalData$GSE65801
#save(CanDataset1, file = "./Data/CanDataSet1.rda")


#CanDataset3 <- getGEOData(c("GSE79973"))
#CanDataset3 <- CanDataset3$originalData$GSE79973

#CanDataset4 <- getGEOData(c("GSE19826"))
#CanDataset4 <- CanDataset4$originalData$GSE19826

#CanDataset5 <- getGEOData(c("GSE49051"))
#CanDataset5 <- CanDataset5$originalData$GSE49051

#save(CanDataset3, file = "./Data/CanDataset3.rda")
#save(CanDataset4, file = "./Data/CanDataset4.rda")
#save(CanDataset5, file = "./Data/CanDataset5.rda")

# save(ValDataset1, ValDataset2, ValDataset3, ValDataset4, file = "./Data/HPyloriCellLineData.rda")

## Load Training data sets (4)
load("./Data/FattyLiverData.rda")

## Load the validation datasets (cell lines)
#load("./Data/HPyloriCellLineData.rda")

## Load the cancer dataset 1
#load("./Data/CanDataSet1.rda")

# load Cancer DataSet2 Pheno
#canPheno2 <- read.delim("./Data/EMTAB1440/EMTAB1440Pheno.txt")

# load Cancer DataSet2 Expr
#canExpr2 <- read.delim("./Data/EMTAB1440/ExpressionNormalized_v2.txt")


#load("./Data/CanDataset3.rda")
#load("./Data/CanDataset4.rda")
#load("./Data/CanDataset5.rda")

###############
#######################################################

## Getting the phenotype data for each data set
pheno1 <- Dataset1$pheno
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno

# Valpheno1 <- ValDataset1$pheno
# Valpheno2 <- ValDataset2$pheno
# Valpheno3 <- ValDataset3$pheno
# Valpheno4 <- ValDataset4$pheno

# canPheno1 <- CanDataset1$pheno 
# 
# canPheno3 <- CanDataset3$pheno
# canPheno4 <- CanDataset4$pheno
# canPheno5 <- CanDataset5$pheno

################################

## Getting the expression data for each data set
## load expr4
#load("/Users/mohamedomar/Documents/Research/Projects/Prostate/Data/expr4.rda")
#ProstateData$originalData$GSE46691$expr <- expr4
expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr

# Valexpr1 <- ValDataset1$expr
# Valexpr2 <- ValDataset2$expr
# Valexpr3 <- ValDataset3$expr
# Valexpr4 <- ValDataset4$expr
# 
# canExpr1 <- CanDataset1$expr
# 
# canExpr3 <- CanDataset3$expr
# canExpr4 <- CanDataset4$expr
# canExpr5 <- CanDataset5$expr

## Checking if the expression data are normalized and log2 transformed
# boxplot(expr1[,1:15], outline= FALSE)
# boxplot(expr2[,1:15], outline= FALSE)
# boxplot(expr3[,1:15], outline = FALSE)
# boxplot(expr4[,1:15], outline= FALSE)
# 
# boxplot(Valexpr1[,1:15], outline = FALSE)
# boxplot(Valexpr2[,1:6], outline = FALSE)
# boxplot(Valexpr3[,1:6], outline = FALSE)
# boxplot(Valexpr4[,1:6], outline = FALSE)
# 
# boxplot(canExpr1[,1:14], outline = FALSE)
#################################################################
## Create a list containing training data sets
AllDataSets <- list(Dataset1, Dataset2, Dataset3, Dataset4)
names(AllDataSets) <- c(Dataset1$formattedName, Dataset2$formattedName, Dataset3$formattedName, Dataset4$formattedName)

###################################################################

## Annotate expression

## Expr1
head(rownames(expr1))
#rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]

#####################
## expr2
# head(rownames(expr2))
# rownames(expr2) <- Dataset2$keys
# expr2 <- expr2[!is.na(rownames(expr2)), ]
# dim(expr2)
# 
# # annotate canExpr2 (illumina)
# tmp <- rownames(expr2)
# expr2$GeneSymbol<- mapIds(org.Hs.eg.db,
#                              keys=tmp,
#                              column="SYMBOL",
#                              keytype="ENTREZID",
#                              multiVals="first")
# 
# canExpr2 <- canExpr2[!duplicated(canExpr2$GeneSymbol), ]
# canExpr2 <- canExpr2[!is.na(canExpr2$GeneSymbol), ]
# 
# rownames(canExpr2) <- canExpr2$GeneSymbol
# canExpr2$Hybridization.REF <- NULL
# canExpr2$GeneSymbol <- NULL
# dim(canExpr2)

#####################
## expr3
head(rownames(expr3))


# annotate expr3
tmp <- expr3$key
expr3$GeneSymbol <- mapIds(org.Hs.eg.db,
                          keys=tmp,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

expr3 <- expr3[!duplicated(expr3$GeneSymbol), ]
expr3 <- expr3[!is.na(expr3$GeneSymbol), ]

rownames(expr3) <- expr3$GeneSymbol
expr3$key <- NULL
expr3$GeneSymbol <- NULL
dim(expr3)

## Filter expr3 by keeping only the genes with cpm > 1 (raw counts > 15) in at least 50% of samples
mycpm <- cpm(expr3)
#plot(expr3_raw[,1],mycpm[,1],xlim=c(0,40),ylim=c(0,3))
#abline(v=15,col=2)
#abline(h=1,col=4)

thresh <- mycpm > 1
keep <- rowSums(thresh) >= ncol(expr3)/2
table(keep)

expr3 <- expr3[keep,]
dim(expr3)

## Visualize density before/after filtering
#plot(density(log2(as.matrix(expr3_raw))))
#plot(density(log2(as.matrix(expr3))))


expr3 <- DGEList(expr3)
#barplot(expr3$samples$lib.size, names.arg = colnames(expr3), las=2)


expr3 <- calcNormFactors(expr3, method = c("TMM"))
#expr3$samples

plotMD(expr3,column=2)
abline(h=0,col="grey")

expr3 <- cpm(expr3, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
#boxplot(expr3[,1:10])

#
# #######################
# expr4
head(rownames(expr4))

# annotate expr4
tmp <- as.character(expr4$entrez_id)
expr4$GeneSymbol <- mapIds(org.Hs.eg.db,
                           keys=tmp,
                           column="SYMBOL",
                           keytype="ENTREZID",
                           multiVals="first")

expr4 <- expr4[!duplicated(expr4$GeneSymbol), ]
expr4 <- expr4[!is.na(expr4$GeneSymbol), ]

rownames(expr4) <- expr4$GeneSymbol
expr4$entrez_id <- NULL
expr4$GeneSymbol <- NULL
dim(expr4)

expr4 <- log2(expr4 + 1)
expr4 <- as.matrix(expr4)
range(expr4)

# #######################
# Valexpr1
# head(rownames(Valexpr1))
# rownames(Valexpr1) <- ValDataset1$keys
# Valexpr1 <- Valexpr1[!is.na(rownames(Valexpr1)), ]
# dim(Valexpr1)
# 
# # #######################
# # Valexpr2
# head(rownames(Valexpr2))
# rownames(Valexpr2) <- ValDataset2$keys
# Valexpr2 <- Valexpr2[!is.na(rownames(Valexpr2)), ]
# dim(Valexpr2)
# 
# # #######################
# # Valexpr3
# head(rownames(Valexpr3))
# rownames(Valexpr3) <- ValDataset3$keys
# Valexpr3 <- Valexpr3[!is.na(rownames(Valexpr3)), ]
# dim(Valexpr3)
# 
# # #######################
# # Valexpr4
# head(rownames(Valexpr4))
# rownames(Valexpr4) <- ValDataset4$keys
# Valexpr4 <- Valexpr4[!is.na(rownames(Valexpr4)), ]
# dim(Valexpr4)
# 
# # #######################
# # cancerExpr1
# head(rownames(canExpr1))
# rownames(canExpr1) <- CanDataset1$keys
# canExpr1 <- canExpr1[!is.na(rownames(canExpr1)), ]
# dim(canExpr1)
# 
# # #######################
# # cancerExpr2
# head(rownames(canExpr2))
# 
# # annotate canExpr2 (illumina)
# tmp <- canExpr2$Hybridization.REF
# canExpr2$GeneSymbol<- mapIds(illuminaHumanv3.db,
#                              keys=tmp,
#                              column="SYMBOL",
#                              keytype="PROBEID",
#                              multiVals="first")
# 
# canExpr2 <- canExpr2[!duplicated(canExpr2$GeneSymbol), ]
# canExpr2 <- canExpr2[!is.na(canExpr2$GeneSymbol), ]
# 
# rownames(canExpr2) <- canExpr2$GeneSymbol
# canExpr2$Hybridization.REF <- NULL
# canExpr2$GeneSymbol <- NULL
# dim(canExpr2)
# 
# # Finally modify the sample names to match those in the phenotype
# colnames(canExpr2) <- gsub("X", "", colnames(canExpr2))
# colnames(canExpr2) <- tolower(colnames(canExpr2))
# colnames(canExpr2)
# 
# # Convert to numeric matrix
# COLS <- colnames(canExpr2)
# ROWS <- rownames(canExpr2)
# canExpr2 <- matrix(as.numeric(unlist(canExpr2)),nrow=nrow(canExpr2))
# rownames(canExpr2) <- ROWS
# colnames(canExpr2) <- COLS
# 
# ###########################
# ## CancerExpr3
# head(rownames(canExpr3))
# rownames(canExpr3) <- CanDataset3$keys
# canExpr3 <- canExpr3[!is.na(rownames(canExpr3)), ]
# dim(canExpr3)
# 
# ###########################
# ## CancerExpr4
# head(rownames(canExpr4))
# rownames(canExpr4) <- CanDataset4$keys
# canExpr4 <- canExpr4[!is.na(rownames(canExpr4)), ]
# dim(canExpr4)
# 
# ###########################
# ## CancerExpr3
# head(rownames(canExpr5))
# rownames(canExpr5) <- CanDataset5$keys
# canExpr5 <- canExpr5[!is.na(rownames(canExpr5)), ]
# dim(canExpr5)

# ############################################################
####################################################################
#### Modify the phenotypes

# Pheno1

pheno1$DiseaseStatus <- pheno1$characteristics_ch1.20
pheno1$DiseaseStatus[pheno1$DiseaseStatus == "nas: Steatosis: 0"] <- "control"
pheno1$DiseaseStatus[pheno1$DiseaseStatus %in% c("nas: Steatosis: 1", "nas: Steatosis: 2", "nas: Steatosis: 3")] <- "case"
pheno1 <- pheno1[-c(62:66), ]
pheno1$DiseaseStatus <- factor(pheno1$DiseaseStatus, levels = c("control","case"))
table(pheno1$DiseaseStatus)

#####
expr1 <- expr1[, colnames(expr1) %in% rownames(pheno1)]
all(rownames(pheno1) == colnames(expr1))
# 
Dataset1$pheno <- pheno1
Dataset1$expr <- expr1
Dataset1$keys <- rownames(expr1)


########## 

## Modify pheno2

# pheno2$DiseaseStatus <- pheno2$`gastritis grade:ch1`
# pheno2$DiseaseStatus[pheno2$DiseaseStatus == "normal"] <- "control"
# pheno2$DiseaseStatus[pheno2$DiseaseStatus %in% c("mild", "IM", "severe")] <- "case"
# pheno2$DiseaseStatus <- factor(pheno2$DiseaseStatus, levels = c("control", "case"))
# table(pheno2$DiseaseStatus)
# 
# # 
# all(rownames(pheno2) == colnames(expr2))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# Dataset2$expr <- expr2
# Dataset2$pheno <- pheno2
# Dataset2$keys <- rownames(expr2)

######
# Modify pheno3
pheno3$DiseaseStatus <- pheno3$`disease:ch1`
pheno3$DiseaseStatus[pheno3$DiseaseStatus %in% c("healthy", "obese")] <- "control"
pheno3$DiseaseStatus[pheno3$DiseaseStatus %in% c("NASH", "NAFLD")] <- "case"

pheno3$DiseaseStatus <- factor(pheno3$DiseaseStatus, levels = c("control", "case"))
table(pheno3$DiseaseStatus)

rownames(pheno3) <- pheno3$description
rownames(pheno3)[32:57] <- paste("X0", rownames(pheno3)[32:57] , sep = "")
rownames(pheno3)[1:31] <- paste("X", rownames(pheno3)[1:31] , sep = "")

pheno3 <- pheno3[order(rownames(pheno3)), ]
expr3 <- expr3[, order(colnames(expr3))]

all(rownames(pheno3) == colnames(expr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset3$pheno <- pheno3
Dataset3$expr <- expr3
Dataset3$keys <- rownames(expr3)

####### 
## Modify pheno4
pheno4$DiseaseStatus <- pheno4$`steatosis grade:ch1`
pheno4$DiseaseStatus[pheno4$DiseaseStatus == 0] <- "control"
pheno4$DiseaseStatus[pheno4$DiseaseStatus %in% c(1,2,3)] <- "case"


pheno4$DiseaseStatus <- factor(pheno4$DiseaseStatus, levels = c("control", "case"))
table(pheno4$DiseaseStatus)

## Modify sample names to match sample names of pheno4
rownames(pheno4) <- pheno4$title
rownames(pheno4) <- paste("X", rownames(pheno4) , sep = "")
all(rownames(pheno4) == colnames(expr4))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset4$expr <- expr4
Dataset4$pheno <- pheno4
Dataset4$keys <- rownames(expr4)

################## 
## Modify valPheno1
# Valpheno1$DiseaseStatus <- Valpheno1$source_name_ch1
# Valpheno1 <- Valpheno1[1:8, ]
# Valpheno1$DiseaseStatus[1:4] <- "case"
# Valpheno1$DiseaseStatus[5:8] <- "control"
# 
# Valpheno1$DiseaseStatus <- factor(Valpheno1$DiseaseStatus, levels = c("control", "case"))
# table(Valpheno1$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno1
# Valexpr1 <- Valexpr1[, colnames(Valexpr1) %in% rownames(Valpheno1)]
# all(rownames(Valpheno1) == colnames(Valexpr1))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# ValDataset1$expr <- Valexpr1
# ValDataset1$pheno <- Valpheno1
# ValDataset1$keys <- rownames(Valexpr1)
# 
# ################## 
# ## Modify valPheno2
# Valpheno2$DiseaseStatus <- Valpheno2$`infection:ch1`
# Valpheno2$DiseaseStatus[Valpheno2$DiseaseStatus == "uninfected"] <- "control"
# Valpheno2$DiseaseStatus[Valpheno2$DiseaseStatus == "H. pylori strain 60190"] <- "case"
# 
# Valpheno2$DiseaseStatus <- factor(Valpheno2$DiseaseStatus, levels = c("control", "case"))
# table(Valpheno2$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno2
# all(rownames(Valpheno2) == colnames(Valexpr2))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# ValDataset2$expr <- Valexpr2
# ValDataset2$pheno <- Valpheno2
# ValDataset2$keys <- rownames(Valexpr2)
# 
# ################## 
# ## Modify valPheno3
# Valpheno3$DiseaseStatus <- Valpheno3$title
# Valpheno3$DiseaseStatus[1:3] <- "control"
# Valpheno3$DiseaseStatus[4:6] <- "case"
# 
# Valpheno3$DiseaseStatus <- factor(Valpheno3$DiseaseStatus, levels = c("control", "case"))
# table(Valpheno3$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno3
# all(rownames(Valpheno3) == colnames(Valexpr3))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# ValDataset3$expr <- Valexpr3
# ValDataset3$pheno <- Valpheno3
# ValDataset3$keys <- rownames(Valexpr3)
# 
# ################## 
# ## Modify valPheno4
# Valpheno4$DiseaseStatus <- Valpheno4$title
# Valpheno4$DiseaseStatus[1:3] <- "control"
# Valpheno4$DiseaseStatus[4:6] <- "case"
# 
# Valpheno4$DiseaseStatus <- factor(Valpheno4$DiseaseStatus, levels = c("control", "case"))
# table(Valpheno4$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno3
# all(rownames(Valpheno4) == colnames(Valexpr4))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# ValDataset4$expr <- Valexpr4
# ValDataset4$pheno <- Valpheno4
# ValDataset4$keys <- rownames(Valexpr4)
# 
# ################## 
# ## Modify cancer pheno1
# canPheno1$DiseaseStatus <- canPheno1$`tissue:ch1`
# canPheno1$DiseaseStatus[canPheno1$DiseaseStatus == "normal gastric tissue"] <- "control"
# canPheno1$DiseaseStatus[canPheno1$DiseaseStatus == "gastric tumor"] <- "case"
# 
# canPheno1$DiseaseStatus <- factor(canPheno1$DiseaseStatus, levels = c("control", "case"))
# table(canPheno1$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno3
# all(rownames(canPheno1) == colnames(canExpr1))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# CanDataset1$expr <- canExpr1
# CanDataset1$pheno <- canPheno1
# CanDataset1$keys <- rownames(canExpr1)
# 
# ################## 
# ## Modify cancer pheno2
# 
# rownames(canPheno2) <- canPheno2$Source.Name
# 
# 
# canPheno2$DiseaseStatus <- canPheno2$Factor.Value.disease.
# canPheno2$DiseaseStatus[canPheno2$DiseaseStatus == "normal"] <- "control"
# canPheno2$DiseaseStatus[canPheno2$DiseaseStatus == "gastric adenocarcinoma"] <- "case"
# 
# canPheno2$DiseaseStatus <- factor(canPheno2$DiseaseStatus, levels = c("control", "case"))
# table(canPheno2$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno3
# all(rownames(canPheno2) == colnames(canExpr2))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# CanDataset2 <- list()
# CanDataset2$expr <- matrix(as.numeric(unlist(canExpr2)),nrow=nrow(canExpr2))
# CanDataset2$pheno <- canPheno2
# CanDataset2$keys <- rownames(canExpr2)
# CanDataset2$formattedName <- "EMTAB1440"
# 
# 
# ################## 
# ## Modify cancer pheno3
# 
# canPheno3$DiseaseStatus <- canPheno3$`tissue:ch1`
# canPheno3$DiseaseStatus[canPheno3$DiseaseStatus == "gastric mucosa"] <- "control"
# canPheno3$DiseaseStatus[canPheno3$DiseaseStatus == "gastric adenocarcinoma"] <- "case"
# 
# canPheno3$DiseaseStatus <- factor(canPheno3$DiseaseStatus, levels = c("control", "case"))
# table(canPheno3$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno3
# all(rownames(canPheno3) == colnames(canExpr3))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# CanDataset3 <- list()
# CanDataset3$expr <- canExpr3
# CanDataset3$pheno <- canPheno3
# CanDataset3$keys <- rownames(canExpr3)
# CanDataset3$formattedName <- "GSE79973"
# 
# ################## 
# ## Modify cancer pheno4
# 
# canPheno4$DiseaseStatus <- canPheno4$`tissue type:ch1`
# canPheno4$DiseaseStatus[canPheno4$DiseaseStatus == "noncancer tissue"] <- "control"
# canPheno4$DiseaseStatus[canPheno4$DiseaseStatus == "normal gastric tissue"] <- "control"
# canPheno4$DiseaseStatus[canPheno4$DiseaseStatus == "gastric cancer tissue"] <- "case"
# 
# canPheno4$DiseaseStatus <- factor(canPheno4$DiseaseStatus, levels = c("control", "case"))
# table(canPheno4$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno4
# all(rownames(canPheno4) == colnames(canExpr4))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# CanDataset4 <- list()
# CanDataset4$expr <- canExpr4
# CanDataset4$pheno <- canPheno4
# CanDataset4$keys <- rownames(canExpr4)
# CanDataset4$formattedName <- "GSE19826"
# 
# ################## 
# ## Modify cancer pheno5
# 
# canPheno5$DiseaseStatus <- canPheno5$`tissue:ch1`
# canPheno5$DiseaseStatus[canPheno5$DiseaseStatus == "normal gastric"] <- "control"
# canPheno5$DiseaseStatus[canPheno5$DiseaseStatus == "gastric cancer tumor"] <- "case"
# 
# canPheno5$DiseaseStatus <- factor(canPheno5$DiseaseStatus, levels = c("control", "case"))
# table(canPheno5$DiseaseStatus)
# 
# ## Modify sample names to match sample names of valpheno5
# all(rownames(canPheno5) == colnames(canExpr5))
# 
# ## Finally, replace the expression and phenotype data in the dataset with the new modified versions
# CanDataset5 <- list()
# CanDataset5$expr <- canExpr5
# CanDataset5$pheno <- canPheno5
# CanDataset5$keys <- rownames(canExpr5)
# CanDataset5$formattedName <- "GSE49051"
# 

#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset1 <- classFunction(Dataset1, column = "DiseaseStatus", diseaseTerms = c("case"))
#Dataset2 <- classFunction(Dataset2, column = "DiseaseStatus", diseaseTerms = c("case"))
Dataset3 <- classFunction(Dataset3, column = "DiseaseStatus", diseaseTerms = c("case"))
Dataset4 <- classFunction(Dataset4, column = "DiseaseStatus", diseaseTerms = c("case"))

# ValDataset1 <- classFunction(ValDataset1, column = "DiseaseStatus", diseaseTerms = c("case"))
# ValDataset2 <- classFunction(ValDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))
# ValDataset3 <- classFunction(ValDataset3, column = "DiseaseStatus", diseaseTerms = c("case"))
# ValDataset4 <- classFunction(ValDataset4, column = "DiseaseStatus", diseaseTerms = c("case"))
# 
# CanDataset1 <- classFunction(CanDataset1, column = "DiseaseStatus", diseaseTerms = c("case"))
# CanDataset2 <- classFunction(CanDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))
# CanDataset3 <- classFunction(CanDataset3, column = "DiseaseStatus", diseaseTerms = c("case"))
# CanDataset4 <- classFunction(CanDataset4, column = "DiseaseStatus", diseaseTerms = c("case"))
# CanDataset5 <- classFunction(CanDataset5, column = "DiseaseStatus", diseaseTerms = c("case"))

############################################################################
############################################################################
#########################################################################
## The metaanalysis

## Creating the meta object
AllDataSets <- list(Dataset1, Dataset3, Dataset4)
names(AllDataSets) <- c(Dataset1$formattedName, Dataset3$formattedName, Dataset4$formattedName)

FattyLiverMeta <- list()
FattyLiverMeta$originalData <- AllDataSets

## Check the meta object before the metaanalysis
checkDataObject(FattyLiverMeta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis

#HPylori_Meta <- geneSymbolCorrection(HPylori_Meta)

## Run the meta analysis
FattyLiver_metaanalysis <- runMetaAnalysis(FattyLiverMeta, runLeaveOneOutAnalysis = F, maxCores = 3)

# effectSize = 0/ FDR = 0.1/ Nstudies = 3

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
FattyLiver_metaanalysis <- filterGenes(FattyLiver_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.05, numberStudiesThresh = 3)

## Assigning a name to the filter
filter <- FattyLiver_metaanalysis$filterResults[[1]]
filter

PositiveGenes <- filter$posGeneNames
PositiveGenes
NegativeGenes <- filter$negGeneNames
NegativeGenes

#save(PositiveGenes, NegativeGenes, file = "./Objs/NewSigGenes_Pre.rda")


## Summarize filter results
filter_summary <- summarizeFilterResults(metaObject = FattyLiver_metaanalysis, getMostRecentFilter(FattyLiver_metaanalysis))

## Save the filter
save(filter, file = "./Objs/filter.rda")

## Save a table of the positive genes and negative genes
#write.table(filter_summary$pos, file = "./Objs/Meta/Positive_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")
#write.table(filter_summary$neg, file = "./Objs/Meta/Negative_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")

## Modify the gene signature for more accuracy and AUC
# Using forward search 
#New_filter <- forwardSearch(metaObject = HPylori_metaanalysis, filterObject = filter)

#save(New_filter, file = "./Objs/NewFilter.rda")
#load("./Objs/NewFilter.rda")

## Replace the old filter with the new smaller one
#HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$posGeneNames <- New_filter$posGeneNames
#HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$negGeneNames <- New_filter$negGeneNames

#New_filter <- HPylori_metaanalysis$filterResults[[1]]
#New_filter_summary <- summarizeFilterResults(metaObject = HPylori_metaanalysis, getMostRecentFilter(HPylori_metaanalysis))

## Save the tables of positive and negative genes
# write.table(New_filter_summary$pos, file = "./Objs/NewFilter_Positive_genes.csv", quote = T, sep = "\t", col.names = T, row.names = T)
# write.table(New_filter_summary$neg, file = "./Objs/NewFilter_Negative_genes.csv", quote = T, sep = "\t", col.names = T, row.names = T)
# 
# ## Gene names
# PositiveGenes <- New_filter$posGeneNames
# PositiveGenes
# NegativeGenes <- New_filter$negGeneNames
# NegativeGenes
# 
# write.table(PositiveGenes, file = "./UpRegulatedGenes.txt")
# write.table(NegativeGenes, file = "./DownRegulatedGenes.txt")

#New_SignatureGenes <- c(PositiveGenes, NegativeGenes)
#save(New_SignatureGenes, file = "./Objs/NewSigGenes.rda")

## Create a summary ROC curve (Training data sets)
set.seed(333)
png(filename = "./Figs/Pooled_ROC_TrainingDatasets.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = FattyLiver_metaanalysis, filterObject = filter)
dev.off()



#### ## Load the new genes from Mohamed
# NewGenes <- read.delim("./Data/FromMohamed/0-Heli_compined_gene_signature.txt")
# ## Replace the old filter with the new smaller one
# HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$posGeneNames <- NewGenes$Upregulated
# HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$negGeneNames <- NewGenes$Downregulated[-c(25:31)]
# 
# New_filter2 <- HPylori_metaanalysis$filterResults[[1]]
# New_filter2_summary <- summarizeFilterResults(metaObject = HPylori_metaanalysis, getMostRecentFilter(HPylori_metaanalysis))
# 
# ## Create a summary ROC curve (Training data sets)
# set.seed(333)
# png(filename = "./Figs/Pooled_ROC_TrainingDatasets_NewFilter2.png", width = 2000, height = 2000, res = 300)
# pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter2)
# dev.off()



##########################################3
library(enrichR)
library(clusterProfiler)

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018", "KEGG_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
Enriched_PositiveGns <- enrichr(genes = PositiveGenes, databases = dbs)
Enriched_NegativeGns <- enrichr(genes = NegativeGenes, databases = dbs)
printEnrich(Enriched_PositiveGns, "PosGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
printEnrich(Enriched_NegativeGns, "NegGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
# 
Pos_GO_BP <- Enriched_PositiveGns["GO_Biological_Process_2018"]
Pos_GO_BP <- Pos_GO_BP$GO_Biological_Process_2018
Pos_GO_BP <- Pos_GO_BP[Pos_GO_BP$P.value <= 0.05, ]
# 
Neg_GO_BP <- Enriched_NegativeGns["GO_Biological_Process_2018"]
Neg_GO_BP <- Neg_GO_BP$GO_Biological_Process_2018
Neg_GO_BP <- Neg_GO_BP[Neg_GO_BP$P.value <= 0.05, ]


Pos_KEGG <- Enriched_PositiveGns["KEGG_2019_Human"]
Pos_KEGG <- Pos_KEGG$KEGG_2019_Human
Pos_KEGG <- Pos_KEGG[Pos_KEGG$P.value <= 0.05, ]

Neg_KEGG <- Enriched_NegativeGns["KEGG_2019_Human"]
Neg_KEGG <- Neg_KEGG$KEGG_2019_Human
Neg_KEGG <- Neg_KEGG[Neg_KEGG$P.value <= 0.05, ]



###########################
Expr1_Up <- t(expr1[PositiveGenes, ])
Data1_Up <- cbind(Expr1_Up, pheno1$DiseaseStatus)
Data1_Up <- as.data.frame(Data1_Up)
colnames(Data1_Up)[63] <- "DiseaseStatus"
Data1_Up$DiseaseStatus <- factor(Data1_Up$DiseaseStatus)
levels(Data1_Up$DiseaseStatus) <- c("control", "case")
table(Data1_Up$DiseaseStatus)
Data1_Up <- Data1_Up[!is.na(Data1_Up$DiseaseStatus), ]

p <- ggplot(Data1_Up, aes(DiseaseStatus, CDKN1A))
p + geom_boxplot(na.rm = T)


#########
expr3_Z <- t(scale(t(expr3), scale = T, center = T))

aPKC_Status <- ifelse(expr3_Z["PRKCI", ] < 0 & expr3_Z["PRKCZ", ] < 0, "low_aPKC", "normal_aPKC")
table(aPKC_Status)
pheno3$aPKC_Status <- aPKC_Status

Expr3_Up <- t(expr3[PositiveGenes, ])
Data3_Up <- as.data.frame(Expr3_Up)
Data3_Up$aPKC_Status <- pheno3$aPKC_Status
Data3_Up$Disease_Status <- pheno3$DiseaseStatus
levels(Data3_Up$Disease_Status) <- c("normal", "NAFLD")
#Data3_Up <- as.data.frame(Data3_Up)
#colnames(Data3_Up)[63] <- "aPKC_Status"
Data3_Up$aPKC_Status <- factor(Data3_Up$aPKC_Status, levels = c("normal_aPKC", "low_aPKC"))
#levels(Data3_Up$aPKC_Status) <- c("control", "case")
table(Data3_Up$aPKC_Status)
Data3_Up <- Data3_Up[!is.na(Data3_Up$aPKC_Status), ]

p <- ggplot(Data3_Up, aes(Disease_Status, BAX, fill = aPKC_Status))
p + geom_boxplot(na.rm = T)


###
Expr3_Down <- t(expr3[NegativeGenes, ])
Data3_Down <- as.data.frame(Expr3_Down)
Data3_Down$aPKC_Status <- pheno3$aPKC_Status
Data3_Down$Disease_Status <- pheno3$DiseaseStatus
levels(Data3_Down$Disease_Status) <- c("normal", "NAFLD")
#Data3_Down <- as.data.frame(Data3_Down)
#colnames(Data3_Down)[63] <- "aPKC_Status"
Data3_Down$aPKC_Status <- factor(Data3_Down$aPKC_Status, levels = c("normal_aPKC", "low_aPKC"))
#levels(Data3_Down$aPKC_Status) <- c("control", "case")
table(Data3_Down$aPKC_Status)
Data3_Down <- Data3_Down[!is.na(Data3_Down$aPKC_Status), ]

p <- ggplot(Data3_Down, aes(Disease_Status, CRADD, fill = aPKC_Status))
p + geom_boxplot(na.rm = T)


#####################
expr4_Z <- t(scale(t(expr4), scale = T, center = T))

aPKC_Status <- ifelse(expr4_Z["PRKCI", ] < 0 & expr4_Z["PRKCZ", ] < 0, "low_aPKC", "normal_aPKC")
table(aPKC_Status)
pheno4$aPKC_Status <- aPKC_Status

Expr4_Up <- t(expr4[PositiveGenes, ])
Data4_Up <- as.data.frame(Expr4_Up)
Data4_Up$aPKC_Status <- pheno4$aPKC_Status
Data4_Up$Disease_Status <- pheno4$DiseaseStatus
levels(Data4_Up$Disease_Status) <- c("normal", "NAFLD")
#Data4_Up <- as.data.frame(Data4_Up)
#colnames(Data4_Up)[64] <- "aPKC_Status"
Data4_Up$aPKC_Status <- factor(Data4_Up$aPKC_Status, levels = c("normal_aPKC", "low_aPKC"))
#levels(Data4_Up$aPKC_Status) <- c("control", "case")
table(Data4_Up$aPKC_Status)
Data4_Up <- Data4_Up[!is.na(Data4_Up$aPKC_Status), ]

p <- ggplot(Data4_Up, aes(Disease_Status, CASP3, fill = aPKC_Status))
p + geom_boxplot(na.rm = T)


###
Expr4_Down <- t(expr4[NegativeGenes, ])
Data4_Down <- as.data.frame(Expr4_Down)
Data4_Down$aPKC_Status <- pheno4$aPKC_Status
Data4_Down$Disease_Status <- pheno4$DiseaseStatus
levels(Data4_Down$Disease_Status) <- c("normal", "NAFLD")
#Data4_Down <- as.data.frame(Data4_Down)
#colnames(Data4_Down)[64] <- "aPKC_Status"
Data4_Down$aPKC_Status <- factor(Data4_Down$aPKC_Status, levels = c("normal_aPKC", "low_aPKC"))
#levels(Data4_Down$aPKC_Status) <- c("control", "case")
table(Data4_Down$aPKC_Status)
Data4_Down <- Data4_Down[!is.na(Data4_Down$aPKC_Status), ]

p <- ggplot(Data4_Down, aes(Disease_Status, CRADD, fill = aPKC_Status))
p + geom_boxplot(na.rm = T)




