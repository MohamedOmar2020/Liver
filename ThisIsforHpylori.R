


## Load the H.Pylori Signature
load("./Objs/HPylori_Signature.rda")


# Validation dataset1
#png(filename = "./Figs/ROC_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
#rocPlot(datasetObject = Dataset1, filterObject = New_filter2)
#dev.off()

# Validation dataset2
pheno3$score <- calculateScore(filterObject = New_filter2, datasetObject = Dataset3)

sscurves_FattyLiver1 <- evalmod(scores = pheno3$score, labels = pheno3$DiseaseStatus)

rocPlot(datasetObject = Dataset3, filterObject = New_filter2) 
FattyLiver1 <- autoplot(sscurves_FattyLiver1, curvetype = c("ROC")) + labs(title = "GSE126848 (Fatty liver disease)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.24"), size =5)

# Validation dataset3
pheno4$score <- calculateScore(filterObject = New_filter2, datasetObject = Dataset4)
sscurves_FattyLiver2 <- evalmod(scores = pheno4$score, labels = pheno4$DiseaseStatus)
rocPlot(datasetObject = Dataset4, filterObject = New_filter2)
FattyLiver2 <- autoplot(sscurves_FattyLiver2, curvetype = c("ROC")) + labs(title = "GSE130970 (Fatty liver disease)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.40"), size =5)


##################################
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

## ROC Plot
PhenoLiver$score <- calculateScore(filterObject = New_filter2, datasetObject = LiverCancer)
sscurves_LiverCancer <- evalmod(scores = PhenoLiver$score, labels = PhenoLiver$DiseaseStatus)
rocPlot(datasetObject = LiverCancer, filterObject = New_filter2)
LiverCancer <- autoplot(sscurves_LiverCancer, curvetype = c("ROC")) + labs(title = "GSE101685 (Hepatocellular carcinoma)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.14"), size =5)

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

## ROC Plot
PhenoLiverAdenoma$score <- calculateScore(filterObject = New_filter2, datasetObject = LiverAdenoma)
sscurves_LiverAdenoma <- evalmod(scores = PhenoLiverAdenoma$score, labels = PhenoLiverAdenoma$DiseaseStatus)
rocPlot(datasetObject = LiverAdenoma, filterObject = New_filter2)
LiverAdenoma <- autoplot(sscurves_LiverAdenoma, curvetype = c("ROC")) + labs(title = "GSE88839 (Liver adenoma)") + annotate("text", x = .75, y = .25, label = paste("AUC = 0.36"), size =5)

#########################################
## Crohn's disease vs normal
CrohnDisease <- getGEOData("GSE83448")
CrohnDisease <- CrohnDisease$originalData$GSE83448


# Modifiy expression
ExprCrohnDisease <- CrohnDisease$expr
head(rownames(ExprCrohnDisease))
rownames(ExprCrohnDisease) <- CrohnDisease$keys
dim(ExprCrohnDisease)
ExprCrohnDisease <- ExprCrohnDisease[!is.na(rownames(ExprCrohnDisease)), ]


# Modify pheno
PhenoCrohnDisease <- CrohnDisease$pheno
PhenoCrohnDisease$DiseaseStatus <- PhenoCrohnDisease$`inflammation:ch1`
PhenoCrohnDisease$DiseaseStatus[PhenoCrohnDisease$DiseaseStatus %in% c("Non-inflamed margin", "Control")] <- "control"
PhenoCrohnDisease$DiseaseStatus[PhenoCrohnDisease$DiseaseStatus == "Inflamed margin"] <- "case"
PhenoCrohnDisease$DiseaseStatus <- factor(PhenoCrohnDisease$DiseaseStatus, levels = c("control","case"))
table(PhenoCrohnDisease$DiseaseStatus)

#####
all(rownames(PhenoCrohnDisease) == colnames(ExprCrohnDisease))
# 
CrohnDisease$pheno <- PhenoCrohnDisease
CrohnDisease$expr <- ExprCrohnDisease
CrohnDisease$keys <- rownames(ExprCrohnDisease)

CrohnDisease <- classFunction(CrohnDisease, column = "DiseaseStatus", diseaseTerms = c("case"))

PhenoCrohnDisease$score <- calculateScore(filterObject = New_filter2, datasetObject = CrohnDisease)
sscurves_Crohns1 <- evalmod(scores = PhenoCrohnDisease$score, labels = PhenoCrohnDisease$DiseaseStatus)
rocPlot(datasetObject = CrohnDisease, filterObject = New_filter2)
Crohns1 <- autoplot(sscurves_Crohns1, curvetype = c("ROC")) + labs(title = "GSE83448 (Crohns disease)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.63"), size =5)

#########################################
## Crohn's disease vs normal
CrohnDisease2 <- getGEOData("GSE112366")
CrohnDisease2 <- CrohnDisease2$originalData$GSE112366

# Modifiy expression
ExprCrohnDisease2 <- CrohnDisease2$expr
head(rownames(ExprCrohnDisease2))
rownames(ExprCrohnDisease2) <- CrohnDisease2$keys
dim(ExprCrohnDisease2)
ExprCrohnDisease2 <- ExprCrohnDisease2[!is.na(rownames(ExprCrohnDisease2)), ]


# Modify pheno
PhenoCrohnDisease2 <- CrohnDisease2$pheno
PhenoCrohnDisease2$DiseaseStatus <- PhenoCrohnDisease2$`diagnosis:ch1`
PhenoCrohnDisease2$DiseaseStatus[PhenoCrohnDisease2$DiseaseStatus == "Normal (non-IBD)"] <- "control"
PhenoCrohnDisease2$DiseaseStatus[PhenoCrohnDisease2$DiseaseStatus == "Crohn's disease (CD)"] <- "case"
PhenoCrohnDisease2$DiseaseStatus <- factor(PhenoCrohnDisease2$DiseaseStatus, levels = c("control","case"))
table(PhenoCrohnDisease2$DiseaseStatus)

#####
all(rownames(PhenoCrohnDisease2) == colnames(ExprCrohnDisease2))
# 
CrohnDisease2$pheno <- PhenoCrohnDisease2
CrohnDisease2$expr <- ExprCrohnDisease2
CrohnDisease2$keys <- rownames(ExprCrohnDisease2)

CrohnDisease2 <- classFunction(CrohnDisease2, column = "DiseaseStatus", diseaseTerms = c("case"))


## ROC Plot
PhenoCrohnDisease2$score <- calculateScore(filterObject = New_filter2, datasetObject = CrohnDisease2)
sscurves_Crohns2 <- evalmod(scores = PhenoCrohnDisease2$score, labels = PhenoCrohnDisease2$DiseaseStatus)
rocPlot(datasetObject = CrohnDisease2, filterObject = New_filter2)
Crohns2 <- autoplot(sscurves_Crohns2, curvetype = c("ROC")) + labs(title = "GSE112366 (Crohns disease)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.665"), size =5)

############################################
#########################################
##  IBD vs normal
# IBD <- getGEOData("GSE75214")
# IBD <- IBD$originalData$GSE75214
# 
# 
# # Modifiy expression
# ExprIBD <- IBD$expr
# head(rownames(ExprIBD))
# rownames(ExprIBD) <- IBD$keys
# dim(ExprIBD)
# ExprIBD <- ExprIBD[!is.na(rownames(ExprIBD)), ]
# 
# 
# # Modify pheno
# PhenoIBD <- IBD$pheno
# PhenoIBD$DiseaseStatus <- PhenoIBD$`disease:ch1`
# PhenoIBD$DiseaseStatus[PhenoIBD$DiseaseStatus == "control"] <- "control"
# PhenoIBD$DiseaseStatus[PhenoIBD$DiseaseStatus %in% c("Crohn's disease", "ulcerative colitis", "CD")] <- "case"
# PhenoIBD$DiseaseStatus <- factor(PhenoIBD$DiseaseStatus, levels = c("control","case"))
# table(PhenoIBD$DiseaseStatus)
# 
# #####
# all(rownames(PhenoIBD) == colnames(ExprIBD))
# # 
# IBD$pheno <- PhenoIBD
# IBD$expr <- ExprIBD
# IBD$keys <- rownames(ExprIBD)
# 
# IBD <- classFunction(IBD, column = "DiseaseStatus", diseaseTerms = c("case"))
# 
# PhenoIBD$score <- calculateScore(filterObject = New_filter2, datasetObject = IBD)
# sscurves_IBD <- evalmod(scores = PhenoIBD$score, labels = PhenoIBD$DiseaseStatus)
# rocPlot(datasetObject = IBD, filterObject = New_filter2)
# 
# #prcPlot(datasetObject = IBD, filterObject = New_filter2)
# 
# #violinPlot(filterObject = New_filter2, datasetObject = IBD, labelColumn = "DiseaseStatus")


########################################
#########################################
##  OA vs normal
OA <- getGEOData("GSE117999")
OA <- OA$originalData$GSE117999


# Modifiy expression
ExprOA <- OA$expr
head(rownames(ExprOA))
rownames(ExprOA) <- OA$keys
dim(ExprOA)
ExprOA <- ExprOA[!is.na(rownames(ExprOA)), ]


# Modify pheno
PhenoOA <- OA$pheno
PhenoOA$DiseaseStatus <- PhenoOA$title
PhenoOA$DiseaseStatus[1:13] <- "control"
PhenoOA$DiseaseStatus[14:24] <- "case"
PhenoOA$DiseaseStatus <- factor(PhenoOA$DiseaseStatus, levels = c("control","case"))
table(PhenoOA$DiseaseStatus)

#####
all(rownames(PhenoOA) == colnames(ExprOA))
# 
OA$pheno <- PhenoOA
OA$expr <- ExprOA
OA$keys <- rownames(ExprOA)

OA <- classFunction(OA, column = "DiseaseStatus", diseaseTerms = c("case"))

## ROC Plot
PhenoOA$score <- calculateScore(filterObject = New_filter2, datasetObject = OA)
sscurves_OA <- evalmod(scores = PhenoOA$score, labels = PhenoOA$DiseaseStatus)
rocPlot(datasetObject = OA, filterObject = New_filter2)
OA <- autoplot(sscurves_OA, curvetype = c("ROC")) + labs(title = "GSE117999 (Osteoarthritis)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.59"), size =5)



###########################################################################3
## Plot the figure


tiff(filename = "./Figs/IntersecSig_OtherInflammatoryDs.tiff", width = 1500, height = 1000)
((FattyLiver1 / FattyLiver2  + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
                  (LiverAdenoma / LiverCancer + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
                  (Crohns1 / Crohns2 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
                  (OA  + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) 
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'Refined Signature performance in other inflammatory diseases',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 17, face = "bold"))
  )
dev.off()


