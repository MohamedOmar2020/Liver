

rm(list = ls())

###############################################################
## Divide the colon cancer dataset according to PRKCI expression

load("./Data/ColoRectalCancerDataset.rda")
load("./Objs/filter.rda")

ValPheno1 <- ColoRectalCancerDataset$pheno

ValExpr1 <- ColoRectalCancerDataset$expr


## ValExpr1
head(rownames(ValExpr1))
rownames(ValExpr1) <- ColoRectalCancerDataset$keys
dim(ValExpr1)
ValExpr1 <- ValExpr1[!is.na(rownames(ValExpr1)), ]

############################################
# Divide the expression matrix into PRKCI low (< 0 ) vs high (> 0)

PRKCI_low <- which(ValExpr1["PRKCI", ] > 0)
PRKCI_low_exp <- ValExpr1[, PRKCI_low]

PRKCI_high_exp <- ValExpr1[, -PRKCI_low]

## Divide the phenotype

#PRKCI Low
PRKCI_low_pheno <- ValPheno1[rownames(ValPheno1) %in% colnames(PRKCI_low_exp), ]

PRKCI_low_pheno$DiseaseStatus <- PRKCI_low_pheno$`tissue type:ch1`
PRKCI_low_pheno$DiseaseStatus[PRKCI_low_pheno$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
PRKCI_low_pheno$DiseaseStatus[PRKCI_low_pheno$DiseaseStatus == "Tumor tissue"] <- "case"
PRKCI_low_pheno$DiseaseStatus <- factor(PRKCI_low_pheno$DiseaseStatus, levels = c("control","case"))
table(PRKCI_low_pheno$DiseaseStatus)


# PRKCI High
PRKCI_high_pheno <- ValPheno1[rownames(ValPheno1) %in% colnames(PRKCI_high_exp), ]

PRKCI_high_pheno$DiseaseStatus <- PRKCI_high_pheno$`tissue type:ch1`
PRKCI_high_pheno$DiseaseStatus[PRKCI_high_pheno$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
PRKCI_high_pheno$DiseaseStatus[PRKCI_high_pheno$DiseaseStatus == "Tumor tissue"] <- "case"
PRKCI_high_pheno$DiseaseStatus <- factor(PRKCI_high_pheno$DiseaseStatus, levels = c("control","case"))
table(PRKCI_high_pheno$DiseaseStatus)

#################################################################
## Make new datasets and test the signature

# PRKCI low
PRKC_low_dataset <- list()
PRKC_low_dataset$pheno <- PRKCI_low_pheno
PRKC_low_dataset$expr <- PRKCI_low_exp
PRKC_low_dataset$keys <- rownames(PRKCI_low_exp)
PRKC_low_dataset$formattedName <- "PRKC_Low"

PRKC_low_dataset <- classFunction(PRKC_low_dataset, column = "DiseaseStatus", diseaseTerms = c("case"))

# Test
rocPlot(datasetObject = PRKC_low_dataset, filterObject = filter)
prcPlot(datasetObject = PRKC_low_dataset, filterObject = filter)
violinPlot(filterObject = filter, datasetObject = PRKC_low_dataset, labelColumn = "DiseaseStatus")

#########################
# PRKCI high
PRKC_high_dataset <- list()
PRKC_high_dataset$pheno <- PRKCI_high_pheno
PRKC_high_dataset$expr <- PRKCI_high_exp
PRKC_high_dataset$keys <- rownames(PRKCI_high_exp)
PRKC_high_dataset$formattedName <- "PRKC_high"

PRKC_high_dataset <- classFunction(PRKC_high_dataset, column = "DiseaseStatus", diseaseTerms = c("case"))

# Test
rocPlot(datasetObject = PRKC_high_dataset, filterObject = filter)
prcPlot(datasetObject = PRKC_high_dataset, filterObject = filter)
violinPlot(filterObject = filter, datasetObject = PRKC_high_dataset, labelColumn = "DiseaseStatus")




