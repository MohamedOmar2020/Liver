
library(dplyr)
library(ArrayTools)

pheno1 <- pheno1[!is.na(pheno1$DiseaseStatus), ]

expr1 <- expr1[,colnames(expr1) %in% rownames(pheno1)]
expr1_GSEA <- as.data.frame(expr1)

expr1_GSEA$NAME <- rownames(expr1_GSEA)

expr1_GSEA$DESCRIBTION <- rep(1, nrow(expr1_GSEA))

expr1_GSEA <- expr1_GSEA %>%
  select(c("NAME", "DESCRIBTION"), everything())



write.table(expr1_GSEA, file="./Text/expr1.txt",row.names = FALSE,sep = "\t", na = "na")

output.cls(pheno1, 78, filename = "./Text/phenotype1")
