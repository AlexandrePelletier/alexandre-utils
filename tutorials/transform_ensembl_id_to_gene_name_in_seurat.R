
#transform ensembl id to hgnc symbol (gene name) in seurat

library(data.table)
library(Seurat)
s1<-readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/simpleaf/sample1.rds')

head(rownames(s1))
# [1] "ENSG00000268020"  "ENSG00000240361"  "ENSG00000186092"  "ENSG00000248527"  "ENSG00000228794" 
# [6] "ENSG00000188157"

#translation betweem ensembl and hgnc symbok
genes_trans<-fread('/projectnb/tcwlab/LabMember/adpelle1/projects/APOE_Jorganoid/outputs/CellRangerCount/sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz',select = 1:2,col.names = c('gene_id','hgnc_symbol'))

#            gene_id hgnc_symbol
# 1: ENSG00000237613     FAM138A
# 2: ENSG00000186092       OR4F5
# 3: ENSG00000238009  AL627309.1
# 4: ENSG00000239945  AL627309.3
# 5: ENSG00000239906  AL627309.2
# ---                            
#   36596: ENSG00000277836  AC141272.1
# 36597: ENSG00000278633  AC023491.2
# 36598: ENSG00000276017  AC007325.1
# 36599: ENSG00000278817  AC007325.4
# 36600: ENSG00000277196  AC007325.2


have.gene_name<-!is.na(genes_trans[rownames(s1),on='gene_id']$hgnc_symbol)
sum(have.gene_name) #34027/36600 have gene names

s1<-s1[have.gene_name,]

mat<-s1[['RNA']]$counts
rownames(mat)<-genes_trans[rownames(mat),on='gene_id']$hgnc_symbol

sum(duplicated(rownames(mat)))#7 duplicated genes

s1[['RNA']]<-CreateAssayObject(counts = mat)
head(rownames(s1)) #"OR4F5"      "LINC01128"  "AGRN"       "MRPL20-AS1" "AL645728.1" "MIB2"   

#keep the ensembl id in metadata one
s1[['RNA']][['gene_id']]<-setNames(genes_trans[rownames(s1),on='hgnc_symbol']$gene_id,
                                                           rownames(s1))
s1[['RNA']]@meta.features
