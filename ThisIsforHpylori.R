


## Load the H.Pylori Signature
load("./Objs/HPylori_Signature.rda")


# Validation dataset1
png(filename = "./Figs/ROC_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = Dataset1, filterObject = New_filter2)
dev.off()

# Validation dataset2

png(filename = "./Figs/ROC_Test2_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = Dataset3, filterObject = New_filter2)
dev.off()

# Validation dataset3

png(filename = "./Figs/ROC_Test3_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = Dataset4, filterObject = New_filter2)
dev.off()


#### Using PRC plot

# val dataset1
png(filename = "./Figs/PRC_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = Dataset1, filterObject = New_filter2)
dev.off()

# val dataset2
png(filename = "./Figs/PRC_Test2_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = Dataset3, filterObject = New_filter2)
dev.off()

# val dataset3
png(filename = "./Figs/PRC_Test3_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = Dataset4, filterObject = New_filter2)
dev.off()

##################################

#### Using violin plot

# Val Dataset1 
png(filename = "./Figs/ViolinPlot_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = Dataset1, labelColumn = "DiseaseStatus")
dev.off()

# Val Dataset2
png(filename = "./Figs/ViolinPlot_Test2_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = Dataset3, labelColumn = "DiseaseStatus")
dev.off()

# Val Dataset3
png(filename = "./Figs/ViolinPlot_Test3_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = Dataset4, labelColumn = "DiseaseStatus")
dev.off()

##############
#### Calculate a signature score (Z score) and add it to the phenotype table

# Val Dataset1
pheno1$score <- calculateScore(filterObject = New_filter2, datasetObject = Dataset1)

# Val Dataset2
pheno3$score <- calculateScore(filterObject = New_filter2, datasetObject = Dataset3)

# Val Dataset3
pheno4$score <- calculateScore(filterObject = New_filter2, datasetObject = Dataset4)


#############################
##### Predictions

# Val Dataset1
Test_Predictions1 <- ifelse(pheno1$score >= -0.5570008, "case", "control")
confusionMatrix(as.factor(Test_Predictions1), pheno1$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions1), actuals = pheno1$DiseaseStatus)

# Val Dataset2
Test_Predictions2 <- ifelse(pheno3$score >= -0.5570008, "case", "control")
confusionMatrix(as.factor(Test_Predictions2), pheno3$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions2), actuals = pheno3$DiseaseStatus)

# Val Dataset3
Test_Predictions3 <- ifelse(pheno4$score >= -0.5570008, "case", "control")
confusionMatrix(as.factor(Test_Predictions3), pheno4$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions3), actuals = pheno4$DiseaseStatus)


###########################################################################3
## Colon cancer
load("./Data/ColoRectalCancerDataset.rda")
#load("./Objs/filter.rda")

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

rocPlot(datasetObject = ColoRectalCancerDataset, filterObject = New_filter2)

prcPlot(datasetObject = ColoRectalCancerDataset, filterObject = New_filter2)


violinPlot(filterObject = New_filter2, datasetObject = ColoRectalCancerDataset, labelColumn = "DiseaseStatus")

# Val Dataset1
ValPheno1$score <- calculateScore(filterObject = New_filter2, datasetObject = ColoRectalCancerDataset)

Test_Predictions1 <- ifelse(ValPheno1$score >= -0.5570008, "case", "control")
confusionMatrix(as.factor(Test_Predictions1), ValPheno1$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions1), actuals = ValPheno1$DiseaseStatus)

###################
## Liver metastasis vs normal liver tissue

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

## Test the signature

#load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = ColoRectalLiverMetsDataset, filterObject = New_filter2)

prcPlot(datasetObject = ColoRectalLiverMetsDataset, filterObject = New_filter2)

violinPlot(filterObject = filter, datasetObject = ColoRectalLiverMetsDataset, labelColumn = "DiseaseStatus")

##################################

## Liver cancer vs normal
LiverCancer <- getGEOData("GSE101685")
LiverCancer <- LiverCancer$originalData$GSE101685


# Modifiy expression
ExprLiver <- LiverCancer$expr
head(rownames(ExprLiver))
rownames(ExprLiver) <- LiverCancer$keys
dim(ExprLiver)
ExprLiver <- ExprLiver[!is.na(rownames(ExprLiver)), ]


# Modify pheno
PhenoLiver <- LiverCancer$pheno
PhenoLiver$DiseaseStatus <- PhenoLiver$`tissue:ch1`
PhenoLiver$DiseaseStatus[PhenoLiver$DiseaseStatus == "Normal liver"] <- "control"
PhenoLiver$DiseaseStatus[PhenoLiver$DiseaseStatus == "Hepatocellular carcinoma"] <- "case"
PhenoLiver$DiseaseStatus <- factor(PhenoLiver$DiseaseStatus, levels = c("control","case"))
table(PhenoLiver$DiseaseStatus)

#####
all(rownames(PhenoLiver) == colnames(ExprLiver))
# 
LiverCancer$pheno <- PhenoLiver
LiverCancer$expr <- ExprLiver
LiverCancer$keys <- rownames(ExprLiver)

LiverCancer <- classFunction(LiverCancer, column = "DiseaseStatus", diseaseTerms = c("case"))

## Test the signature

#load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = LiverCancer, filterObject = New_filter2)

prcPlot(datasetObject = LiverCancer, filterObject = New_filter2)

violinPlot(filterObject = New_filter2, datasetObject = LiverCancer, labelColumn = "DiseaseStatus")

##################################

## Liver adenoma vs normal
LiverAdenoma <- getGEOData("GSE88839")
LiverAdenoma <- LiverAdenoma$originalData$GSE88839


# Modifiy expression
ExprLiverAdenoma <- LiverAdenoma$expr
head(rownames(ExprLiverAdenoma))
rownames(ExprLiverAdenoma) <- LiverAdenoma$keys
dim(ExprLiverAdenoma)
ExprLiverAdenoma <- ExprLiverAdenoma[!is.na(rownames(ExprLiverAdenoma)), ]


# Modify pheno
PhenoLiverAdenoma <- LiverAdenoma$pheno
PhenoLiverAdenoma$DiseaseStatus <- PhenoLiverAdenoma$`disease state:ch1`
PhenoLiverAdenoma$DiseaseStatus[PhenoLiverAdenoma$DiseaseStatus == "Non tumor liver"] <- "control"
PhenoLiverAdenoma$DiseaseStatus[PhenoLiverAdenoma$DiseaseStatus == "Solid Tumor"] <- "case"
PhenoLiverAdenoma$DiseaseStatus <- factor(PhenoLiverAdenoma$DiseaseStatus, levels = c("control","case"))
table(PhenoLiverAdenoma$DiseaseStatus)

#####
all(rownames(PhenoLiverAdenoma) == colnames(ExprLiverAdenoma))
# 
LiverAdenoma$pheno <- PhenoLiverAdenoma
LiverAdenoma$expr <- ExprLiverAdenoma
LiverAdenoma$keys <- rownames(ExprLiverAdenoma)

LiverAdenoma <- classFunction(LiverAdenoma, column = "DiseaseStatus", diseaseTerms = c("case"))

## Test the signature

#load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = LiverAdenoma, filterObject = New_filter2)

prcPlot(datasetObject = LiverAdenoma, filterObject = New_filter2)

violinPlot(filterObject = New_filter2, datasetObject = LiverAdenoma, labelColumn = "DiseaseStatus")





