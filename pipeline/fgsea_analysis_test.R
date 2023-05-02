
#installation
require("devtools")
devtools::install_github("ctlab/fgsea")



#analysis
library(fgsea)
library(data.table)



pathways_dir<-'/projectnb/tcwlab/MSigDB/' #folder containing all  reference pathways files in gmt format can be obtained at: https://drive.google.com/drive/folders/1-9CxYMLg7UBhzwigAIN8ACd8sDnHxUzE?usp=sharing

diffexp_file_path<-"/projectnb/tcwlab/LabMember/adpelle1/projects/juao_autophagy/outputs/1-APOE4_astro_human/res_pseudobulkDESeq2_Dementia_APOE4_vs_3_astro_Cov_braak_stage_sex_PMI_age_education_avg.mt.csv.gz"

pathways_to_test<-c("GO_all", "GO_bp", "GO_cc", "GO_mf", 'CP_all',"CP_biocarta", "CP_kegg", "CP_reactome", "CP_pid") # choose pathway gene sets reference you want to test or let like that test for all

GenesInPathways<-function(pathways,source,path_to_pathways_dir=pathways_dir){
  #return list of genes for each pathway
  gmt_mtd<-fread(file.path(path_to_pathways_dir,'gmt_metadata.csv')) 
  
  return(gmtPathways(file.path(path_to_pathways_dir,  gmt_mtd[name==source]$gmt))[pathways])
}

#calculate your gene stats
res_diffexp<-fread(diffexp_file_path)

res_diffexp[,gene_stat:=-log10(pvalue)*sign(log2FoldChange)] #change 'avg_log2FC' for your Fold change column name and 'p_val' for your pvalue column name

#here the gene stat choose to rank genes for fgsea is the -log10(pvalue) multiply by 1 or -1 depending of the sens of the change (upregulate or downregulate in your tested condition, respectively)

gene_stats<-res_diffexp[order(gene_stat)]$gene_stat #get the vectors of gene statistics
names(gene_stats)<-res_diffexp[order(gene_stat)]$hgnc_symbol #name these vectors with the gene name (in hgnc symbol format)

#load pathways metadata
gmt_mtd<-fread(file.path(pathways_dir,'gmt_metadata.csv')) 

#run fgsea for every pathway source
res_gsea<-Reduce(rbind,lapply(pathways_to_test, function(p){
  message(paste('testing enrichement for',gmt_mtd[name==p]$desc))
  pathways<- gmtPathways(file.path(pathways_dir,gmt_mtd[name==p]$gmt))
  
  res<-fgsea(pathways,
             stats=gene_stats,minSize=10,maxSize=2000,scoreType='std',eps=0)
  
  return(res[,source:=p])
  
}))


res_gsea[padj<0.05] #display all pathways enriched in up or downregulated genes at adjusted pvalue 0.05
res_gsea[padj<0.05&source=='GO_bp'] #display GO Biological Process enriched in up or downregulated genes at adjusted pvalue 0.05

res_gsea[padj<0.05&NES>0] #display all pathways enriched in upregulated genes 
res_gsea[padj<0.05&NES>0] #display GO Biological Process enriched in upregulated genes 

res_gsea[padj<0.05&NES<0] #display all GO terms enriched in downregulated genes 
#...

#save results
fwrite(res_gsea,"res_fgsea.csv")

#one putative visualisation 
pathways_to_visualize<-'GO_bp' #copy and change to make the plot it for all pathways
topPathwaysUp <- res_gsea[source==pathways_to_visualize][ES > 0][head(order(pval), n=10), pathway] 
topPathwaysDown <- res_gsea[source==pathways_to_visualize][ES < 0][head(order(pval), n=10), pathway] 

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GenesInPathways(topPathways,pathways_to_visualize), gene_stats, res_gsea[source==pathways_to_visualize], gseaParam = 0.5, colwidths = c(20, 6, 2, 2, 2))





