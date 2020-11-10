
## Liver metastasis vs normal liver tissue

rm(list = ls())
#ColoRectalLiverMetsDataset <- getGEOData("GSE38174")
#ColoRectalLiverMetsDataset <- ColoRectalLiverMetsDataset$originalData$GSE38174

#save(ColoRectalLiverMetsDataset, file = "./Data/ColoRectalLiverMetsDataset.rda")
load("./Data/ColoRectalLiverMetsDataset.rda")

ValPheno1 <- ColoRectalLiverMetsDataset$pheno

ValExpr1 <- ColoRectalLiverMetsDataset$expr


## ValExpr1
head(rownames(ValExpr1))
rownames(ValExpr1) <- ColoRectalLiverMetsDataset$keys
dim(ValExpr1)
ValExpr1 <- ValExpr1[!is.na(rownames(ValExpr1)), ]


# ValPheno1

ValPheno1$DiseaseStatus <- ValPheno1$`tumour stage:ch1`
ValPheno1$DiseaseStatus[ValPheno1$DiseaseStatus == "None"] <- "control"
ValPheno1$DiseaseStatus[ValPheno1$DiseaseStatus %in% c("Metasis", "Halo of Metasis")] <- "case"
ValPheno1$DiseaseStatus <- factor(ValPheno1$DiseaseStatus, levels = c("control","case"))
table(ValPheno1$DiseaseStatus)

#####
all(rownames(ValPheno1) == colnames(ValExpr1))
# 
ColoRectalLiverMetsDataset$pheno <- ValPheno1
ColoRectalLiverMetsDataset$expr <- ValExpr1
ColoRectalLiverMetsDataset$keys <- rownames(ValExpr1)

ColoRectalLiverMetsDataset <- classFunction(ColoRectalLiverMetsDataset, column = "DiseaseStatus", diseaseTerms = c("case"))

#################################################
## Test the signature

load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = ColoRectalLiverMetsDataset, filterObject = filter)

prcPlot(datasetObject = ColoRectalLiverMetsDataset, filterObject = filter)

violinPlot(filterObject = filter, datasetObject = ColoRectalLiverMetsDataset, labelColumn = "DiseaseStatus")
