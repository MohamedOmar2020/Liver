
## colon cancer vs normal colon tissue 

rm(list = ls())

#ColoRectalCancerDataset <- getGEOData("GSE74602")
#ColoRectalCancerDataset <- ColoRectalCancerDataset$originalData$GSE74602
#save(ColoRectalCancerDataset, file = "./Data/ColoRectalCancerDataset.rda")

load("./Data/ColoRectalCancerDataset.rda")
load("./Objs/filter.rda")

ValPheno1 <- ColoRectalCancerDataset$pheno

ValExpr1 <- ColoRectalCancerDataset$expr


## ValExpr1
head(rownames(ValExpr1))
rownames(ValExpr1) <- ColoRectalCancerDataset$keys
dim(ValExpr1)
ValExpr1 <- ValExpr1[!is.na(rownames(ValExpr1)), ]


# ValPheno1

ValPheno1$DiseaseStatus <- ValPheno1$`tissue type:ch1`
ValPheno1$DiseaseStatus[ValPheno1$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
ValPheno1$DiseaseStatus[ValPheno1$DiseaseStatus == "Tumor tissue"] <- "case"
ValPheno1$DiseaseStatus <- factor(ValPheno1$DiseaseStatus, levels = c("control","case"))
table(ValPheno1$DiseaseStatus)

#####
all(rownames(ValPheno1) == colnames(ValExpr1))
# 
ColoRectalCancerDataset$pheno <- ValPheno1
ColoRectalCancerDataset$expr <- ValExpr1
ColoRectalCancerDataset$keys <- rownames(ValExpr1)

ColoRectalCancerDataset <- classFunction(ColoRectalCancerDataset, column = "DiseaseStatus", diseaseTerms = c("case"))

rocPlot(datasetObject = ColoRectalCancerDataset, filterObject = filter)

prcPlot(datasetObject = ColoRectalCancerDataset, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = ColoRectalCancerDataset, labelColumn = "DiseaseStatus")
