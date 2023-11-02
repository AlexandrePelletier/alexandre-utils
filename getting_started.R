
library(data.table)
source('visualisation.R')

#Network visualisation of fgsea results using emmaplot function ####
#load fgsea results
res_fgsea<-fread('examples_data/res_fgsea_apoe22_vs_33_microglia.csv.gz')

#filter the datatable to show only the top pathways enriched
#top80
res_fgsea_top<-head(res_fgsea[order(pval)],80)

emmaplot(res_fgsea_top)

#you can reduce the pathway label size
emmaplot(res_fgsea_top,label.size = 2)


#to display more variety of pathways (coming from different source),the results can be prefiltered for
#top pathways by pathways source (subcategory)
res_gsea_topbycat<-res_fgsea_top[padj<0.05, #select only pathways those with padj<0.05
                                 .SD[head(order(padj),30)], #filter the data.table by source subcategory by keeping only the top30 by subcat
                                 by='subcat']

emmaplot(head(res_gsea_topbycat[order(pval)],80), #subset the top80 pathways according to the pval of the enrichment
         label.size = 2)

#this is a ggplot format so the plot can be saved using
ggsave('examples_data/emmaplot_top80_fgsea_apoe22_vs_33_microglia.pdf')


#[New] You can decided to display only the pathway of a certain gene (This gene is in the leading edge genes of the pathway)
emmaplot(res_fgsea_top,show_pathway_of = 'CDK2')



