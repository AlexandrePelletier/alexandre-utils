

#installation
require("devtools")
devtools::install_github("ctlab/fgsea")

#analysis
library(fgsea)
library(data.table)

gos<-fread("goterms.csv") #file available here: https://www.dropbox.com/s/st7k7e3415sb6t7/goterms.csv?dl=0 


res_diffexp<-fread("your_differential_expression_results_file.csv")

res_diffexp[,gene_stat:=-log10(p_val)*sign(avg_log2FC)] #change 'avg_log2FC' for your Fold change column name and 'p_val' for your pvalue column name

#here the score choose to rank genes for fgsea is the -log10(pvalue) multiply by 1 or -1 depending of the sens of the change (upregulate or downregulate in your tested condition, respectively)

gene_stats<-res_diffexp[order(gene_stat)]$gene_stat #get the vectors of gene statistics
names(gene_stats)<-res_diffexp[order(gene_stat)]$gene #name these vectors with the gene name

#perform the fgsea analysis for each GO terms category (MF, BP, CC)
res_gsea<-Reduce(rbind,lapply(split(gos,by="source"),function(go){
  res<-fgsea(split(go$gene,go$description),
             stats=gene_stats,minSize=10,maxSize=2000,scoreType='std',nPermSimple = 10000)
  
  return(res[,source:=unique(go$source)])
  
}))

res_gsea<-merge(res_gsea[,description:=pathway],unique(gos[,.(term,term.size,description)])) #add description of go terms

res_gsea[padj<0.05] #display all GO terms enriched in up or downregulated genes at adjusted pvalue 0.05

res_gsea[padj<0.05&NES>0] #display all GO terms enriched in upregulated genes 

res_gsea[padj<0.05&NES<0] #display all GO terms enriched in downregulated genes 

#visualisation 
gos_list<-split(gos$gene,gos$description)
gos_list

topPathwaysUp <- res_gsea[ES > 0][head(order(pval), n=10), description]
topPathwaysDown <- res_gsea[ES < 0][head(order(pval), n=10), description]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gos_list[topPathways], gene_stats, res_gsea, 
              gseaParam=0.5)

#save results
fwrite(res_gsea,"res_fgsea.csv")
