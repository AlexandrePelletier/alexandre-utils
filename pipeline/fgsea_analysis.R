
#installation
require("devtools")
devtools::install_github("ctlab/fgsea")


#load packages
library(fgsea)
library(data.table)

pathways_dir<-'/projectnb2/tcwlab/MSigDB/' #folder containing all  reference pathways files in gmt format can be obtained at: https://drive.google.com/drive/folders/1-9CxYMLg7UBhzwigAIN8ACd8sDnHxUzE?usp=sharing

diffexp_file_path<-"path/to/your_differential_expression_results_file.csv"

category_to_test<-c("GO", "CP")# choose pathway gene sets reference you want to test or let like that to test for all

#some useful function
GenesInPathways<-function(pathways,cat,pathways_dir){  #return the list of genes for each pathway

  gmt_mtd<-fread(file.path(pathways_dir,'gmt_metadata2.csv')) #need to have the gmt metadata file in the folder containing the pathways
  
  return(gmtPathways(file.path(path_to_pathways_dir,  gmt_mtd[category==cat]$gmt_file))[pathways])
}


#analysis
#calculate your gene stats
res_diffexp<-fread(diffexp_file_path)

res_diffexp[,gene_stat:=-log10(p_val)*sign(avg_log2FC)] #change 'avg_log2FC' for your Fold change column name and 'p_val' for your pvalue column name

#here the gene stat choose to rank genes for fgsea is the -log10(pvalue) multiply by 1 or -1 depending of the sens of the change (upregulate or downregulate in your tested condition, respectively)

gene_stats<-res_diffexp[order(gene_stat)]$gene_stat #get the vectors of gene statistics
names(gene_stats)<-res_diffexp[order(gene_stat)]$gene #name these vectors with the gene name (in hgnc symbol format)

#load pathways metadata
gmt_mtd<-fread(file.path(pathways_dir,'gmt_metadata.csv')) 


#run fgsea for every pathway source

res_gsea<-Reduce(rbind,lapply(category_to_test, function(cat){
  
  message(paste('testing enrichement for',cat,'category'))
  
  pathways<- gmtPathways(file.path(pathways_dir,gmt_mtd[category==cat]$gmt_file))
  
  res<-fgsea(pathways,
             stats=gene_stats,minSize=10,maxSize=2000,scoreType='std',eps=0)
  
  return(res[,category:=cat])
  
}))

un
res_gsea[padj<0.05] #display all pathways enriched in up or downregulated genes at adjusted pvalue 0.05
res_gsea[padj<0.05&source=='GO'] #display GO tenrs enriched in up or downregulated genes at adjusted pvalue 0.05

res_gsea[padj<0.05&NES>0] #display all pathways enriched in upregulated genes 
res_gsea[padj<0.05&NES>0&source=='CP'] #display Cannonical Pathways enriched in upregulated genes 
res_gsea[padj<0.05&NES<0&source=='GO'] #display GO terms enriched in downregulated genes 

#...
#add subcategory pathway size info
pathways_info<-fread(file.path(pathways_dir,'all_CPandGOs_genesets_metadata.csv.gz'))
res_gsea<-merge(res_gsea,pathways_info,by=c('category','pathway'))[order(source,pval)]



#save results


fwrite(res_gsea,"res_fgsea.csv")

#one putative visualisation 
category_to_visualize<-'GO' #copy and change to make the plot it for all pathways
topPathwaysUp <- res_gsea[category==category_to_visualize][ES > 0][head(order(pval), n=10), pathway] 
topPathwaysDown <- res_gsea[source==category_to_visualize][ES < 0][head(order(pval), n=10), pathway] 

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GenesInPathways(topPathways,pathways_to_visualize), gene_stats, res_gsea[source==pathways_to_visualize], gseaParam = 0.5, colwidths = c(20, 6, 2, 2, 2))




