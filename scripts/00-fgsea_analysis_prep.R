library(data.table)
library(fgsea)

#[edit 05/02/23' new references gmt files]


cps<-fgsea::gmtPathways('/projectnb/tcwlab/MSigDB/c2.cp.v2023.1.Hs.symbols.gmt')
pathway_names<-names(cps)
#get the different pathways subsets (derived from different source)
table(sapply(pathway_names, function(x)strsplit(x,'_')[[1]][1]))
#classfied by subcategory: take only the one with >5 pathways, the others are classified as 'others'
subcat_pathways<-names(table(sapply(pathway_names, function(x)strsplit(x,'_')[[1]][1])))[table(sapply(pathway_names, function(x)strsplit(x,'_')[[1]][1]))>5]
cps_mtd<-data.table(
  gmt_file='c2.cp.v2023.1.Hs.symbols.gmt',
  subcat=c(subcat_pathways,'OTHERS'),
  category='CP'
  
)

#same for GO
gos<-fgsea::gmtPathways('/projectnb/tcwlab/MSigDB/c5.go.v2023.1.Hs.symbols.gmt')
go_names<-names(gos)
#get the different go category subsets 
table(sapply(go_names, function(x)strsplit(x,'_')[[1]][1]))

#classfied by subcategory
subcat_gos<-names(table(sapply(go_names, function(x)strsplit(x,'_')[[1]][1])))
gos_mtd<-data.table(
  gmt_file='c5.go.v2023.1.Hs.symbols.gmt',
  subcat=subcat_gos,
  category='GO'
  
)
#merge and save 
gmt_mtd<-rbind(cps_mtd,gos_mtd)
fwrite(gmt_mtd,'/projectnb/tcwlab/MSigDB/gmt_metadata2.csv')

#by pathways/gene set metadata with gene set size
genesets_mtd<-rbind(do.call(rbind,lapply(names(cps),function(p)data.table(pathway=p,gene=cps[[p]],pathway.size=length(cps[[p]]),category='CP',subcat=ifelse(str_extract(p,'^[A-Z]+')%in%gmt_mtd$subcat,str_extract(p,'^[A-Z]+'),'OTHERS')))),
                    do.call(rbind,lapply(names(gos),function(p)data.table(pathway=p,gene=gos[[p]],pathway.size=length(gos[[p]]),category='GO',subcat=ifelse(str_extract(p,'^[A-Z]+')%in%gmt_mtd$subcat,str_extract(p,'^[A-Z]+'),'OTHERS')))))


unique(genesets_mtd[,.(category,subcat,pathway,pathway.size)])

fwrite(unique(genesets_mtd[,.(category,subcat,pathway,pathway.size)]),'/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')


#add mouse gene and gene_id
msig<-genesets_mtd
mouse<-unique(fread('/projectnb/tcwlab/RefData/gencode/HOM_MouseHumanSequence_translation_reformat.csv.gz'),by='Symbol.human')
msig[,mouse_gene:=mouse[gene,on='Symbol.human']$Symbol.mouse]

ensembl<-fread('/projectnb/tcwlab/RefData/gencode/hg38/gencode.v45.gene_id_names.tsv')|>unique(by='gene_name')

msig[,gene_id:=ensembl[gene,on='gene_name']$gene_id]
fwrite(msig,'/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

