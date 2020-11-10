

###
TCGA_Colon_expr2 <- read.delim("./Data/TCGA_Colon2/data_RNA_Seq_v2_expression_RSEM_UQ_Log2.txt")
TCGA_Colon_Pheno2 <- read.delim("./Data/TCGA_Colon2/Clinical.tsv")

#######
## Expression
TCGA_Colon_expr2 <- TCGA_Colon_expr2[!duplicated(TCGA_Colon_expr2$Hugo_Symbol), ]
TCGA_Colon_expr2 <- TCGA_Colon_expr2[!is.na(TCGA_Colon_expr2$Hugo_Symbol), ]
rownames(TCGA_Colon_expr2) <- TCGA_Colon_expr2$Hugo_Symbol
TCGA_Colon_expr2$Hugo_Symbol <- NULL
TCGA_Colon_expr2$Entrez_Gene_Id <- NULL
TCGA_Colon_expr2 <- TCGA_Colon_expr2[!is.na(rownames(TCGA_Colon_expr2)), ]
TCGA_Colon_expr2 <- TCGA_Colon_expr2[!(rownames(TCGA_Colon_expr2) == ""), ]

sel <- which(apply(TCGA_Colon_expr2, 1, function(x) all(is.finite(x)) ))
TCGA_Colon_expr2 <- TCGA_Colon_expr2[sel, ]

#TCGA_Colon_expr2 <- log2(TCGA_Colon_expr2+1)
#C <- as.matrix(expr7)
#plot(density(log2(C)), xlim=c(-5,20))


################
# Pheno

## EMT and hypermutated VS Epithelial
rownames(TCGA_Colon_Pheno2) <- TCGA_Colon_Pheno2$Patient.ID
rownames(TCGA_Colon_Pheno2) <- paste("X", rownames(TCGA_Colon_Pheno2), sep = "")

TCGA_Colon_Pheno2$DiseaseStatus <- TCGA_Colon_Pheno2$INTEGRATED_PHENOTYPE
TCGA_Colon_Pheno2$DiseaseStatus[TCGA_Colon_Pheno2$DiseaseStatus == "Epithelial"] <- "control"
TCGA_Colon_Pheno2$DiseaseStatus[TCGA_Colon_Pheno2$DiseaseStatus %in% c("EMT", "Hypermutated")] <- "case"
TCGA_Colon_Pheno2 <- TCGA_Colon_Pheno2[!is.na(TCGA_Colon_Pheno2$DiseaseStatus), ]
TCGA_Colon_Pheno2$DiseaseStatus <- factor(TCGA_Colon_Pheno2$DiseaseStatus, levels = c("control","case"))
table(TCGA_Colon_Pheno2$DiseaseStatus)

####
## MSI-High vs MSS
TCGA_Colon_Pheno2$DiseaseStatus2 <- TCGA_Colon_Pheno2$MSI_STATUS
TCGA_Colon_Pheno2$DiseaseStatus2[TCGA_Colon_Pheno2$DiseaseStatus2 == "MSS"] <- "control"
TCGA_Colon_Pheno2$DiseaseStatus2[TCGA_Colon_Pheno2$DiseaseStatus2  == "MSI-H"] <- "case"
TCGA_Colon_Pheno2 <- TCGA_Colon_Pheno2[!is.na(TCGA_Colon_Pheno2$DiseaseStatus2), ]
TCGA_Colon_Pheno2$DiseaseStatus2 <- factor(TCGA_Colon_Pheno2$DiseaseStatus2, levels = c("control","case"))
table(TCGA_Colon_Pheno2$DiseaseStatus2)

#####
TCGA_Colon_expr2 <- TCGA_Colon_expr2[, colnames(TCGA_Colon_expr2) %in% rownames(TCGA_Colon_Pheno2)]
TCGA_Colon_Pheno2 <- TCGA_Colon_Pheno2[rownames(TCGA_Colon_Pheno2) %in% colnames(TCGA_Colon_expr2), ]

all(rownames(TCGA_Colon_Pheno2) == colnames(TCGA_Colon_expr2))
# 
TCGA_Colon_expr2 <- as.matrix(TCGA_Colon_expr2)
TCGA_ColonCanDataset2 <- list()

TCGA_ColonCanDataset2$pheno <- TCGA_Colon_Pheno2
TCGA_ColonCanDataset2$expr <- TCGA_Colon_expr2
TCGA_ColonCanDataset2$keys <- rownames(TCGA_Colon_expr2)
TCGA_ColonCanDataset2$formattedName <- "TCGA_Colon"

TCGA_ColonCanDataset2 <- classFunction(TCGA_ColonCanDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))

rocPlot(datasetObject = TCGA_ColonCanDataset2, filterObject = filter)

prcPlot(datasetObject = TCGA_ColonCanDataset2, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = TCGA_ColonCanDataset2, labelColumn = "DiseaseStatus")

############3
TCGA_ColonCanDataset2 <- classFunction(TCGA_ColonCanDataset2, column = "DiseaseStatus2", diseaseTerms = c("case"))

rocPlot(datasetObject = TCGA_ColonCanDataset2, filterObject = filter)

prcPlot(datasetObject = TCGA_ColonCanDataset2, filterObject = filter)


violinPlot(filterObject = filter, datasetObject = TCGA_ColonCanDataset2, labelColumn = "DiseaseStatus2")

#######################################################
TCGA_Colon_expr2 <- t(scale(t(TCGA_Colon_expr2), scale = T, center = T))

PRKCI_low <- which(TCGA_Colon_expr2["PRKCI", ] < 0 & TCGA_Colon_expr2["PRKCZ", ] < 0)
PRKCI_low_exp <- TCGA_Colon_expr2[, PRKCI_low]

PRKCI_high_exp <- TCGA_Colon_expr2[, -PRKCI_low]

## Divide the phenotype

#PRKCI Low
PRKCI_low_pheno <- TCGA_Colon_Pheno2[rownames(TCGA_Colon_Pheno2) %in% colnames(PRKCI_low_exp), ]

#PRKCI_low_pheno$DiseaseStatus <- PRKCI_low_pheno$DiseaseStatus`
#PRKCI_low_pheno$DiseaseStatus[PRKCI_low_pheno$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
#PRKCI_low_pheno$DiseaseStatus[PRKCI_low_pheno$DiseaseStatus == "Tumor tissue"] <- "case"
#PRKCI_low_pheno$DiseaseStatus <- factor(PRKCI_low_pheno$DiseaseStatus, levels = c("control","case"))
table(PRKCI_low_pheno$DiseaseStatus)
table(PRKCI_low_pheno$DiseaseStatus2)


# PRKCI High
PRKCI_high_pheno <- TCGA_Colon_Pheno2[rownames(TCGA_Colon_Pheno2) %in% colnames(PRKCI_high_exp), ]

#PRKCI_high_pheno$DiseaseStatus <- PRKCI_high_pheno$`tissue type:ch1`
#PRKCI_high_pheno$DiseaseStatus[PRKCI_high_pheno$DiseaseStatus == "Normal non-cancerous tissue"] <- "control"
#PRKCI_high_pheno$DiseaseStatus[PRKCI_high_pheno$DiseaseStatus == "Tumor tissue"] <- "case"
#PRKCI_high_pheno$DiseaseStatus <- factor(PRKCI_high_pheno$DiseaseStatus, levels = c("control","case"))
table(PRKCI_high_pheno$DiseaseStatus)
table(PRKCI_high_pheno$DiseaseStatus2)

#################################################################
## Make new datasets and test the signature

# PRKCI low
PRKC_low_dataset <- list()
PRKC_low_dataset$pheno <- PRKCI_low_pheno
PRKC_low_dataset$expr <- PRKCI_low_exp
PRKC_low_dataset$keys <- rownames(PRKCI_low_exp)
PRKC_low_dataset$formattedName <- "PRKC_Low"

PRKC_low_dataset <- classFunction(PRKC_low_dataset, column = "DiseaseStatus", diseaseTerms = c("case"))
PRKC_low_dataset <- classFunction(PRKC_low_dataset, column = "DiseaseStatus2", diseaseTerms = c("case"))

# Test
rocPlot(datasetObject = PRKC_low_dataset, filterObject = filter)
prcPlot(datasetObject = PRKC_low_dataset, filterObject = filter)
violinPlot(filterObject = filter, datasetObject = PRKC_low_dataset, labelColumn = "DiseaseStatus2")

#########################
# PRKCI high
PRKC_high_dataset <- list()
PRKC_high_dataset$pheno <- PRKCI_high_pheno
PRKC_high_dataset$expr <- PRKCI_high_exp
PRKC_high_dataset$keys <- rownames(PRKCI_high_exp)
PRKC_high_dataset$formattedName <- "PRKC_high"

PRKC_high_dataset <- classFunction(PRKC_high_dataset, column = "DiseaseStatus", diseaseTerms = c("case"))
PRKC_high_dataset <- classFunction(PRKC_high_dataset, column = "DiseaseStatus2", diseaseTerms = c("case"))

# Test
rocPlot(datasetObject = PRKC_high_dataset, filterObject = filter)
prcPlot(datasetObject = PRKC_high_dataset, filterObject = filter)
violinPlot(filterObject = filter, datasetObject = PRKC_high_dataset, labelColumn = "DiseaseStatus")



################
#TCGA_Colon_expr2_Z <- t(scale(t(TCGA_Colon_expr2), scale = T, center = T))

aPKC_Status <- ifelse(TCGA_Colon_expr2["PRKCI", ] < 0 & TCGA_Colon_expr2["PRKCZ", ] < 0, "low_aPKC", "normal_aPKC")
table(aPKC_Status)
TCGA_Colon_Pheno2$aPKC_Status <- aPKC_Status

TCGA_Colon_expr2_up <- t(TCGA_Colon_expr2[PositiveGenes, ])
Data_Up <- as.data.frame(TCGA_Colon_expr2_up)
Data_Up$aPKC_Status <- TCGA_Colon_Pheno2$aPKC_Status
Data_Up$Disease_Status <- TCGA_Colon_Pheno2$DiseaseStatus
levels(Data_Up$Disease_Status) <- c("MSS", "MSI-H")
#Data3_Up <- as.data.frame(Data3_Up)
#colnames(Data3_Up)[63] <- "aPKC_Status"
Data_Up$aPKC_Status <- factor(Data_Up$aPKC_Status, levels = c("normal_aPKC", "low_aPKC"))
#levels(Data3_Up$aPKC_Status) <- c("control", "case")
table(Data_Up$aPKC_Status)
Data_Up <- Data_Up[!is.na(Data_Up$aPKC_Status), ]

p <- ggplot(Data_Up, aes(Disease_Status, IL7, fill = aPKC_Status))
p + geom_boxplot(na.rm = T)


###
TCGA_Colon_expr2_Down <- t(TCGA_Colon_expr2[NegativeGenes %in% rownames(TCGA_Colon_expr2), ])
Data_Down <- as.data.frame(TCGA_Colon_expr2_Down)
Data_Down$aPKC_Status <- TCGA_Colon_Pheno2$aPKC_Status
Data_Down$Disease_Status <- TCGA_Colon_Pheno2$DiseaseStatus
levels(Data_Down$Disease_Status) <- c("MSS", "MSI-High")
#Data3_Down <- as.data.frame(Data3_Down)
#colnames(Data3_Down)[63] <- "aPKC_Status"
Data_Down$aPKC_Status <- factor(Data_Down$aPKC_Status, levels = c("normal_aPKC", "low_aPKC"))
#levels(Data3_Down$aPKC_Status) <- c("control", "case")
table(Data_Down$aPKC_Status)
Data_Down <- Data_Down[!is.na(Data_Down$aPKC_Status), ]

p <- ggplot(Data_Down, aes(Disease_Status, ICAM3, fill = aPKC_Status))
p + geom_boxplot(na.rm = T)


