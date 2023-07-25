
library(data.table)
source('visualisation.R')

#Network visualisation of fgsea results using emmaplot function####
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
                                 .SD[head(order(padj),30)],
                                 by='subcat']

emmaplot(head(res_gsea_topbycat[order(pval)],80),label.size = 2)
#this is a ggplot format so the plot can be saved using
ggsave('examples_data/emmaplot_top80_fgsea_apoe22_vs_33_microglia.pdf')

