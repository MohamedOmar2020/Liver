
## Liver metastasis vs normal liver tissue 2

rm(list = ls())

# ColoRectalLiverMetsDataset2 <- getGEOData("GSE40367")
# ColoRectalLiverMetsDataset2 <- ColoRectalLiverMetsDataset2$originalData$GSE40367
# save(ColoRectalLiverMetsDataset2, file = "./Data/ColoRectalLiverMetsDataset2.rda")

load("./Data/ColoRectalLiverMetsDataset2.rda")
load("./Objs/filter.rda")

ValPheno2 <- ColoRectalLiverMetsDataset2$pheno

ValExpr2 <- ColoRectalLiverMetsDataset2$expr


## ValExpr2
head(rownames(ValExpr2))
rownames(ValExpr2) <- ColoRectalLiverMetsDataset2$keys
dim(ValExpr2)
ValExpr2 <- ValExpr2[!is.na(rownames(ValExpr2)), ]


# ValPheno2

ValPheno2$DiseaseStatus <- ValPheno2$`tissue:ch1`
ValPheno2$DiseaseStatus[ValPheno2$DiseaseStatus == "tumor endothelium from normal liver"] <- "control"
ValPheno2$DiseaseStatus[ValPheno2$DiseaseStatus == "tumor endothelium from liver metastasis"] <- "case"
ValPheno2 <- ValPheno2[ValPheno2$DiseaseStatus %in% c("control", "case"), ]
ValPheno2$DiseaseStatus <- factor(ValPheno2$DiseaseStatus, levels = c("control","case"))
table(ValPheno2$DiseaseStatus)

#####
ValExpr2 <- ValExpr2[, colnames(ValExpr2) %in% rownames(ValPheno2)]
all(rownames(ValPheno2) == colnames(ValExpr2))
# 
ColoRectalLiverMetsDataset2$pheno <- ValPheno2
ColoRectalLiverMetsDataset2$expr <- ValExpr2
ColoRectalLiverMetsDataset2$keys <- rownames(ValExpr2)

ColoRectalLiverMetsDataset2 <- classFunction(ColoRectalLiverMetsDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))

rocPlot(datasetObject = ColoRectalLiverMetsDataset2, filterObject = filter)

prcPlot(datasetObject = ColoRectalLiverMetsDataset2, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = ColoRectalLiverMetsDataset2, labelColumn = "DiseaseStatus")

######################################

## Using both the primary tumor and metastasis as "case" VS normal liver tissue as "control
#ColonCancerDataset <- getGEOData("GSE40367")
#ColonCancerDataset <- ColonCancerDataset$originalData$GSE40367

load("./Data/ColoRectalLiverMetsDataset2.rda")


ValPheno2 <- ColoRectalLiverMetsDataset2$pheno

ValExpr2 <- ColoRectalLiverMetsDataset2$expr


## ValExpr2
head(rownames(ValExpr2))
rownames(ValExpr2) <- ColoRectalLiverMetsDataset2$keys
dim(ValExpr2)
ValExpr2 <- ValExpr2[!is.na(rownames(ValExpr2)), ]


# ValPheno2

ValPheno2$DiseaseStatus <- ValPheno2$`tissue:ch1`
ValPheno2$DiseaseStatus[ValPheno2$DiseaseStatus == "tumor endothelium from normal liver"] <- "control"
ValPheno2$DiseaseStatus[ValPheno2$DiseaseStatus %in% c("tumor endothelium from liver metastasis", "tumor endothelium from colon adenocarcinoma tissue")] <- "case"
ValPheno2 <- ValPheno2[ValPheno2$DiseaseStatus %in% c("control", "case"), ]
ValPheno2$DiseaseStatus <- factor(ValPheno2$DiseaseStatus, levels = c("control","case"))
table(ValPheno2$DiseaseStatus)

#####
ValExpr2 <- ValExpr2[, colnames(ValExpr2) %in% rownames(ValPheno2)]
all(rownames(ValPheno2) == colnames(ValExpr2))
# 
ColoRectalLiverMetsDataset2$pheno <- ValPheno2
ColoRectalLiverMetsDataset2$expr <- ValExpr2
ColoRectalLiverMetsDataset2$keys <- rownames(ValExpr2)

ColoRectalLiverMetsDataset2 <- classFunction(ColoRectalLiverMetsDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))

rocPlot(datasetObject = ColoRectalLiverMetsDataset2, filterObject = filter)

prcPlot(datasetObject = ColoRectalLiverMetsDataset2, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = ColoRectalLiverMetsDataset2, labelColumn = "DiseaseStatus")

###################################################
# Comparing HCC to normal liver

#ColonCancerDataset <- getGEOData("GSE40367")
#ColonCancerDataset <- ColonCancerDataset$originalData$GSE40367

load("./Data/ColoRectalLiverMetsDataset2.rda")


ValPheno2 <- ColoRectalLiverMetsDataset2$pheno

ValExpr2 <- ColoRectalLiverMetsDataset2$expr


## ValExpr2
head(rownames(ValExpr2))
rownames(ValExpr2) <- ColoRectalLiverMetsDataset2$keys
dim(ValExpr2)
ValExpr2 <- ValExpr2[!is.na(rownames(ValExpr2)), ]


# ValPheno2

ValPheno2$DiseaseStatus <- ValPheno2$`tissue:ch1`
ValPheno2$DiseaseStatus[ValPheno2$DiseaseStatus == "tumor endothelium from normal liver"] <- "control"
ValPheno2$DiseaseStatus[ValPheno2$DiseaseStatus == "tumor endothelium from HCC tumor tissue"] <- "case"
ValPheno2 <- ValPheno2[ValPheno2$DiseaseStatus %in% c("control", "case"), ]
ValPheno2$DiseaseStatus <- factor(ValPheno2$DiseaseStatus, levels = c("control","case"))
table(ValPheno2$DiseaseStatus)

#####
ValExpr2 <- ValExpr2[, colnames(ValExpr2) %in% rownames(ValPheno2)]
all(rownames(ValPheno2) == colnames(ValExpr2))
# 
ColoRectalLiverMetsDataset2$pheno <- ValPheno2
ColoRectalLiverMetsDataset2$expr <- ValExpr2
ColoRectalLiverMetsDataset2$keys <- rownames(ValExpr2)

ColoRectalLiverMetsDataset2 <- classFunction(ColoRectalLiverMetsDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))

rocPlot(datasetObject = ColoRectalLiverMetsDataset2, filterObject = filter)

prcPlot(datasetObject = ColoRectalLiverMetsDataset2, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = ColoRectalLiverMetsDataset2, labelColumn = "DiseaseStatus")
