
## High vs low fat diet
library("readxl")

rm(list = ls())

#ColoRectalLiverMetsDataset3 <- getGEOData("GSE151165")
#ColoRectalLiverMetsDataset3 <- ColoRectalLiverMetsDataset3$originalData$GSE151165

#save(ColoRectalLiverMetsDataset3, file = "./Data/ColoRectalLiverMetsDataset3.rda")

load("./Data/ColoRectalLiverMetsDataset3.rda")

ValPheno3 <- ColoRectalLiverMetsDataset3$pheno

ValExpr3 <- read_excel("./Data/GSE151165_RNA_seq.raw_read_count.xlsx")

ValExpr3 <- ValExpr3[,-32]
## ValExpr1
head(rownames(ValExpr3))
Genes <- ValExpr3$NAME
rownames(ValExpr3) <- ValExpr3$NAME
ValExpr3$NAME <- NULL

# Convert to numeric matrix
COLS <- colnames(ValExpr3)
ROWS <- Genes
ValExpr3 <- matrix(as.numeric(unlist(ValExpr3)),nrow=nrow(ValExpr3))
rownames(ValExpr3) <- ROWS
colnames(ValExpr3) <- COLS

thresh <- mycpm > 1
keep <- rowSums(thresh) >= ncol(ValExpr3)/2
table(keep)

ValExpr3 <- ValExpr3[keep,]
dim(ValExpr3)

ValExpr3 <- DGEList(ValExpr3)

ValExpr3 <- calcNormFactors(ValExpr3, method = c("TMM"))

plotMD(ValExpr3,column=2)
abline(h=0,col="grey")

ValExpr3 <- cpm(ValExpr3, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)


# ValPheno1

ValPheno3$DiseaseStatus <- ValPheno3$`tissue:ch1`
ValPheno3$DiseaseStatus[ValPheno3$DiseaseStatus == "Normal"] <- "control"
ValPheno3$DiseaseStatus[ValPheno3$DiseaseStatus == "Tumor"] <- "case"
ValPheno3$DiseaseStatus <- factor(ValPheno3$DiseaseStatus, levels = c("control","case"))
table(ValPheno3$DiseaseStatus)

#####
rownames(ValPheno3) <- ValPheno3$title
rownames(ValPheno3) <- gsub("\\_.+", "", rownames(ValPheno3))

ValPheno3 <- ValPheno3[order(rownames(ValPheno3)), ]
ValExpr3 <- ValExpr3[, order(colnames(ValExpr3))]

all(rownames(ValPheno3) == colnames(ValExpr3))

# 
ColoRectalLiverMetsDataset3$pheno <- ValPheno3
ColoRectalLiverMetsDataset3$expr <- ValExpr3
ColoRectalLiverMetsDataset3$keys <- rownames(ValExpr3)

ColoRectalLiverMetsDataset3 <- classFunction(ColoRectalLiverMetsDataset3, column = "DiseaseStatus", diseaseTerms = c("case"))

#################################################
## Test the signature

load("./Objs/filter.rda")

## ROC Plot
rocPlot(datasetObject = ColoRectalLiverMetsDataset3, filterObject = filter)

prcPlot(datasetObject = ColoRectalLiverMetsDataset3, filterObject = filter)

violinPlot(filterObject = filter, datasetObject = ColoRectalLiverMetsDataset3, labelColumn = "DiseaseStatus")
