

TCGA_Colon_expr <- read.delim("./Data/TCGA_Colon/data_RNA_Seq_v2_expression_median.txt")
TCGA_Colon_Pheno <- read.delim("./Data/TCGA_Colon/Clinical.tsv")

#######
## Expression
TCGA_Colon_expr <- TCGA_Colon_expr[!duplicated(TCGA_Colon_expr$Hugo_Symbol), ]
TCGA_Colon_expr <- TCGA_Colon_expr[!is.na(TCGA_Colon_expr$Hugo_Symbol), ]
rownames(TCGA_Colon_expr) <- TCGA_Colon_expr$Hugo_Symbol
TCGA_Colon_expr$Hugo_Symbol <- NULL
TCGA_Colon_expr$Entrez_Gene_Id <- NULL
TCGA_Colon_expr <- TCGA_Colon_expr[!is.na(rownames(TCGA_Colon_expr)), ]
TCGA_Colon_expr <- TCGA_Colon_expr[!(rownames(TCGA_Colon_expr) == ""), ]

sel <- which(apply(TCGA_Colon_expr, 1, function(x) all(is.finite(x)) ))
TCGA_Colon_expr <- TCGA_Colon_expr[sel, ]
TCGA_Colon_expr <- log2(TCGA_Colon_expr+1)
#C <- as.matrix(expr7)
#plot(density(log2(C)), xlim=c(-5,20))


################
# Pheno
rownames(TCGA_Colon_Pheno) <- TCGA_Colon_Pheno$Sample.ID
rownames(TCGA_Colon_Pheno) <- gsub("\\-", "\\.", rownames(TCGA_Colon_Pheno))
TCGA_Colon_Pheno$DiseaseStatus <- TCGA_Colon_Pheno$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code
TCGA_Colon_Pheno$DiseaseStatus[TCGA_Colon_Pheno$DiseaseStatus %in% c("M0", "MX")] <- "control"
TCGA_Colon_Pheno$DiseaseStatus[TCGA_Colon_Pheno$DiseaseStatus %in% c("M1", "M1A", "M1B")] <- "case"
TCGA_Colon_Pheno <- TCGA_Colon_Pheno[!is.na(TCGA_Colon_Pheno$DiseaseStatus), ]
TCGA_Colon_Pheno$DiseaseStatus <- factor(TCGA_Colon_Pheno$DiseaseStatus, levels = c("control","case"))
table(TCGA_Colon_Pheno$DiseaseStatus)

#####
TCGA_Colon_expr <- TCGA_Colon_expr[, colnames(TCGA_Colon_expr) %in% rownames(TCGA_Colon_Pheno)]
TCGA_Colon_Pheno <- TCGA_Colon_Pheno[rownames(TCGA_Colon_Pheno) %in% colnames(TCGA_Colon_expr), ]

TCGA_Colon_expr <- TCGA_Colon_expr[, order(colnames(TCGA_Colon_expr))]
TCGA_Colon_Pheno <- TCGA_Colon_Pheno[order(rownames(TCGA_Colon_Pheno)), ]

all(rownames(TCGA_Colon_Pheno) == colnames(TCGA_Colon_expr))
# 
TCGA_Colon_expr <- as.matrix(TCGA_Colon_expr)
TCGA_ColonCanDataset <- list()

TCGA_ColonCanDataset$pheno <- TCGA_Colon_Pheno
TCGA_ColonCanDataset$expr <- TCGA_Colon_expr
TCGA_ColonCanDataset$keys <- rownames(TCGA_Colon_expr)
TCGA_ColonCanDataset$formattedName <- "TCGA_Colon"

TCGA_ColonCanDataset <- classFunction(TCGA_ColonCanDataset, column = "DiseaseStatus", diseaseTerms = c("case"))

rocPlot(datasetObject = TCGA_ColonCanDataset, filterObject = filter)

prcPlot(datasetObject = TCGA_ColonCanDataset, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = TCGA_ColonCanDataset, labelColumn = "DiseaseStatus")


#######################################################
TCGA_Colon_expr <- t(scale(t(TCGA_Colon_expr), scale = T, center = T))

PRKCI_low <- which(TCGA_Colon_expr["PRKCI", ] < 0 & TCGA_Colon_expr["PRKCI", ] < 0)
PRKCI_low_exp <- TCGA_Colon_expr[, PRKCI_low]

PRKCI_high_exp <- TCGA_Colon_expr[, -PRKCI_low]

## Divide the phenotype

#PRKCI Low
PRKCI_low_pheno <- TCGA_Colon_Pheno[rownames(TCGA_Colon_Pheno) %in% colnames(PRKCI_low_exp), ]

PRKCI_low_pheno$DiseaseStatus <- PRKCI_low_pheno$DiseaseStatus`
#PRKCI_low_pheno$DiseaseStatus[PRKCI_low_pheno$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
#PRKCI_low_pheno$DiseaseStatus[PRKCI_low_pheno$DiseaseStatus == "Tumor tissue"] <- "case"
#PRKCI_low_pheno$DiseaseStatus <- factor(PRKCI_low_pheno$DiseaseStatus, levels = c("control","case"))
table(PRKCI_low_pheno$DiseaseStatus)


# PRKCI High
PRKCI_high_pheno <- TCGA_Colon_Pheno[rownames(TCGA_Colon_Pheno) %in% colnames(PRKCI_high_exp), ]

#PRKCI_high_pheno$DiseaseStatus <- PRKCI_high_pheno$`tissue type:ch1`
#PRKCI_high_pheno$DiseaseStatus[PRKCI_high_pheno$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
#PRKCI_high_pheno$DiseaseStatus[PRKCI_high_pheno$DiseaseStatus == "Tumor tissue"] <- "case"
#PRKCI_high_pheno$DiseaseStatus <- factor(PRKCI_high_pheno$DiseaseStatus, levels = c("control","case"))
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





