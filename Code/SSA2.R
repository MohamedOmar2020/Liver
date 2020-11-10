## SSA vs normal colon and conventional adenomas


rm(list = ls())

SSA_Dataset2 <- getGEOData("GSE117606")
SSA_Dataset2 <- SSA_Dataset2$originalData$GSE117606

save(SSA_Dataset2, file = "./Data/SSA_Dataset2.rda")

load("./Data/SSA_Dataset2.rda")

SSA_pheno2 <- SSA_Dataset2$pheno

SSA_expr2 <- SSA_Dataset2$expr

#SSA_expr <- SSA_expr[,-32]

## SSA_expr
head(rownames(SSA_expr2))
rownames(SSA_expr2) <- SSA_Dataset2$keys
dim(SSA_expr2)
SSA_expr2 <- SSA_expr2[!is.na(rownames(SSA_expr2)), ]


# ValPheno2

SSA_pheno2$DiseaseStatus <- SSA_pheno2$source_name_ch2
SSA_pheno2$DiseaseStatus[SSA_pheno2$DiseaseStatus %in% c("pooled normal colonic tissue", "Conventional polyp")] <- "control"
SSA_pheno2$DiseaseStatus[SSA_pheno2$DiseaseStatus == "Sessile serrated adenoma"] <- "case"
SSA_pheno2$DiseaseStatus <- factor(SSA_pheno2$DiseaseStatus, levels = c("control","case"))
table(SSA_pheno2$DiseaseStatus)

#####

all(rownames(SSA_pheno2) == colnames(SSA_expr2))

# 
SSA_Dataset2$pheno <- SSA_pheno2
SSA_Dataset2$expr <- SSA_expr2
SSA_Dataset2$keys <- rownames(SSA_expr2)

SSA_Dataset2 <- classFunction(SSA_Dataset2, column = "DiseaseStatus", diseaseTerms = c("case"))

#################################################
## Test the signature

load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = SSA_Dataset2, filterObject = filter)

prcPlot(datasetObject = SSA_Dataset2, filterObject = filter)

violinPlot(filterObject = filter, datasetObject = SSA_Dataset2, labelColumn = "DiseaseStatus")
