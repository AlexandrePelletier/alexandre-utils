#melt the matrix expression and metatadata to save a single csv file by gene containing samples metadata
out<-'outputs/02-csv_files_tmm_genes'
dir.create(out)

library(data.table)
source('r_utils.R')

expr<-fread('data/tmm_data.csv')

mtd<-fread('data/tmm_meta.csv')

#melt the matrix (wide to long transformation)
expr_dt<-melt(expr,id.vars = 'gname',variable.name = 'id',value.name = 'value')

#merge with metadata
expr_dt<-merge(expr_dt,
               mtd[,.(id,SEX,cell_type,APOE,indiv_id,proportion)]) #keep only useful metadata

#split bu genes and save
expr_dt_list<-split(expr_dt,by='gname')

lapply(names(expr_dt_list), function(gene_name){
  fwrite(expr_dt_list[[gene_name]],file.path(out,paste0('tmm_',gene_name,'.csv')))
})


#check some gene
expr_dt[,genotype:=ps('APOE',str_extract(id,'44|33|KO'))]
expr_dt[,celltype:=ifelse(str_detect(id,'^LA'),'Ast',ifelse(str_detect(id,'^N'),'Neu',"Mic"))]
Lcn2 
unique(expr_dt$gname)
ggplot(expr_dt[gname%in%c('SERPINA3','LCN2','NGAL')][str_detect(id,'LA_TCW1')][genotype!='APOEKO'])+geom_col(aes(x=id,
                                                                                                          y=value,fill=genotype))+
  facet_wrap('gname')


#use raw data
library(DESeq2)
#population
mtd_rna<-fread('/projectnb/tcwlab/Project/PopulationAPOE/JuliaNote/Alignment_final/featureCounts/RNAseqPhen_aug2021.csv')


expr_rna<-as.matrix(data.frame(fread('/projectnb/tcwlab/Project/PopulationAPOE/JuliaNote/Alignment_final/featureCounts/featureCounts_clean.txt'),row.names = 'V1'))

dds <- DESeqDataSetFromMatrix(expr_rna, data.frame(mtd_rna,row.names = 'sample_id'), design = ~RIN+SEX+disease+A_P_O_E)  
dds <- estimateSizeFactors(dds)

norm_count<- counts(dds, normalized=TRUE)

#save
fwrite(data.table(norm_count,keep.rownames = 'gene_id'),fp(out,'DESEQ2_normalized_count_RNA_population.csv.gz'))

expr_rna_dt<-melt(data.table(norm_count,keep.rownames = 'gene_id'),id.vars = 'gene_id',variable.name = 'id',value.name = 'norm_expr')

expr_rna_dta<-merge(expr_rna_dt,mtd_rna)
expr_rna_dta<-merge(expr_rna_dta[,gene_id:=str_remove(gene_id,'\\.[0-9]+')],fread('~/tcwlab/RefData/gencode/hg38/gencode.v45.gene_id_names.tsv'),by='gene_id')

ggplot(expr_rna_dta[gene_name%in%c('IRF1')][A_P_O_E!='K_O'][,id:=str_remove(id,'[LA]+_')])+
  geom_col(aes(x=id,y=norm_expr,fill=A_P_O_E))+facet_wrap('CellType',scales = 'free_x')+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))+ggtitle('IRF1')



mtdind<-fread('../ipsc_metadata/metadata_population_rna_seq.csv')
mtdind[id=='N7_33']

expr_rna_dta<-merge(expr_rna_dta,mtdind[,.(id,ind_id)],by='id')
ggplot(expr_rna_dta[gene_name%in%c('IRF1')][A_P_O_E!='K_O'][,id:=str_remove(id,'[LA]+_')])+
  geom_col(aes(x=ind_id,y=norm_expr,fill=A_P_O_E))+facet_wrap('CellType',scales = 'free_x')+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))+ggtitle('IRF1')



ggplot(expr_rna_dta[gene_name%in%c('IRF1')][A_P_O_E!='K_O'][CellType=='Neuron'][,id:=str_remove(id,'[LA]+_')])+
  geom_col(aes(x=ind_id,y=norm_expr,fill=A_P_O_E))+facet_wrap('CellType',scales = 'free_x')+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))+
  ggtitle('IRF1')

cpm<-edgeR::cpm(expr_rna)[,mtd_rna[CellType=='Neuron']$id]
cpm[str_extract(rownames(cpm),'ENSG[0-9]+')=='ENSG00000125347']
#isogenic#isoge../ipsc_metadata/metadata_population_rna_seq.csvnic

expr_rna<-as.matrix(data.frame(fread('/projectnb/tcwlab/Project/IsogenicAPOE/Alignment/featureCounts/featureCounts_clean_v3.txt'),row.names = 'V1'))
dim(expr_rna)
expr_rna[1:10,1:10]

mtd_rna<-fread('/projectnb/tcwlab/Project/IsogenicAPOE/Alignment/featureCounts/RNAseqPhen_iso_v3.csv')



dds <- DESeqDataSetFromMatrix(expr_rna, data.frame(mtd_rna,row.names = 'id'), design = ~RIN+A_P_O_E)  
dds <- estimateSizeFactors(dds)

norm_count<- counts(dds, normalized=TRUE)


#save
fwrite(data.table(norm_count,keep.rownames = 'gene_id'),fp(out,'DESEQ2_normalized_count_RNA_isogenics.csv.gz'))




expr_rna_dt<-melt(data.table(norm_count,keep.rownames = 'gene_id'),id.vars = 'gene_id',variable.name = 'id',value.name = 'norm_expr')

expr_rna_dta<-merge(expr_rna_dt,mtd_rna)
expr_rna_dta<-merge(expr_rna_dta[,gene_id:=str_remove(gene_id,'\\.[0-9]+')],fread('~/tcwlab/RefData/gencode/hg38/gencode.v45.gene_id_names.tsv'),by='gene_id')

ggplot(expr_rna_dta[gene_name%in%c('SERPINA3','LCN2','NGAL')][A_P_O_E!='K_O'][CellType=='Astrocyte'][,id:=str_remove(id,'[LA]+_')])+
  geom_col(aes(x=id,y=log1p(norm_expr),fill=A_P_O_E))+facet_grid(gene_name~individual,scales = 'free_x')+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))

