
# Getting started with some utils functions 
## emmaplot function 

For Network visualisation of fgsea results using emmaplot function   

First, load some fgsea results
```r
res_fgsea<-fread('examples_data/res_fgsea_apoe22_vs_33_microglia.csv.gz')
```
This data.table should have at least the columns `pathway`, `NES` and `leadingEdge`.  

filter it to show only the top pathways that you want to display.
For example the top80 enriched based on pvalue
```r
res_fgsea_top<-head(res_fgsea[order(pval)],80)
```
Plot these pathways with emmaplot
```r
emmaplot(res_fgsea_top)
```

you can reduce the pathway label size
```r
emmaplot(res_fgsea_top,label.size = 2)
```
### display more diverse pathways
to display more variety of pathways (coming from different source),the results can be prefiltered for top pathways by pathways source (subcategory)
```r
res_gsea_topbycat<-res_fgsea_top[padj<0.05, #select only pathways those with padj<0.05
                                 .SD[head(order(padj),30)], #filter the data.table by source subcategory by keeping only the top30 by subcat
                                 by='subcat']
```

```
emmaplot(head(res_gsea_topbycat[order(pval)],80), #subset the top80 pathways according to the pval of the enrichment
         label.size = 2)
```
this is a ggplot format so the plot can be saved using
```r
ggsave('examples_data/emmaplot_top80_fgsea_apoe22_vs_33_microglia.pdf')

```
You can see it in `examples_data` folder. here is how it look like:
![emmaplot](examples_data/emmaplot_top80_fgsea_apoe22_vs_33_microglia.pdf)
### Filter pathways by a certain gene
You can decided to display only the pathway of a certain gene (i.e this gene is in the leading edge genes of the pathway)
```r
emmaplot(res_fgsea_top,show_pathway_of = 'CDK2')
```
### Get Only the pathways links
For other purpose, you extract the pathways links by using the GetPathwaysLinks function
```r
pathways_links<-GetPathwaysLinks(res_fgsea_top)
pathways_links

```

|                                                                                pathway1  |                                      pathway2  |
| ---------------------------------------------------------------------------------------- | ---------------------------------------------- |
|  1: GOMF_RNA_POLYMERASE_II_TRANSCRIPTION_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING | GOMF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY |
|  2:                                                                   GOCC_MITOCHONDRION |          GOBP_SMALL_MOLECULE_METABOLIC_PROCESS |
|  3:                           REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY |     REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION |
| 4:                           REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY  |            WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS   |

