## SSA vs normal colon and conventional adenomas
 
library("readxl")

rm(list = ls())

#SSA_Dataset1 <- getGEOData("GSE12514")
#SSA_Dataset1 <- SSA_Dataset1$originalData$GSE12514

#save(SSA_Dataset1, file = "./Data/SSA_Dataset1.rda")

load("./Data/SSA_Dataset1.rda")

SSA_pheno1 <- SSA_Dataset1$pheno

SSA_expr1 <- SSA_Dataset1$expr

#SSA_expr <- SSA_expr[,-32]

## SSA_expr
head(rownames(SSA_expr1))
rownames(SSA_expr1) <- SSA_Dataset1$keys
dim(SSA_expr1)
SSA_expr1 <- SSA_expr1[!is.na(rownames(SSA_expr1)), ]


# ValPheno1

SSA_pheno1$DiseaseStatus <- SSA_pheno1$source_name_ch1
SSA_pheno1$DiseaseStatus[SSA_pheno1$DiseaseStatus %in% c("pooled normal colonic tissue", "Conventional polyp")] <- "control"
SSA_pheno1$DiseaseStatus[SSA_pheno1$DiseaseStatus == "Sessile serrated adenoma"] <- "case"
SSA_pheno1$DiseaseStatus <- factor(SSA_pheno1$DiseaseStatus, levels = c("control","case"))
table(SSA_pheno1$DiseaseStatus)

#####

all(rownames(SSA_pheno1) == colnames(SSA_expr1))

# 
SSA_Dataset1$pheno <- SSA_pheno1
SSA_Dataset1$expr <- SSA_expr1
SSA_Dataset1$keys <- rownames(SSA_expr1)

SSA_Dataset1 <- classFunction(SSA_Dataset1, column = "DiseaseStatus", diseaseTerms = c("case"))

#################################################
## Test the signature

load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = SSA_Dataset1, filterObject = filter)

prcPlot(datasetObject = SSA_Dataset1, filterObject = filter)

violinPlot(filterObject = filter, datasetObject = SSA_Dataset1, labelColumn = "DiseaseStatus")
