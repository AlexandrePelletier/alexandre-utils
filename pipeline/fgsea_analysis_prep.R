library(data.table)
library(fgsea)

gmt_mtd <- data.table(
  name=c("GO_all", "GO_bp", "GO_cc", "GO_mf", "CP_biocarta", "CP_kegg", "CP_reactome", "CP_pid"),
  
  
  gmt=c("c5.go.v2023.1.Hs.symbols.gmt"
        ,"c5.go.bp.v2023.1.Hs.symbols.gmt" 
        ,"c5.go.cc.v2023.1.Hs.symbols.gmt" 
        ,"c5.go.mf.v2023.1.Hs.symbols.gmt"
        ,"c2.cp.biocarta.v2023.1.Hs.symbols.gmt"
        ,"c2.cp.kegg.v2023.1.Hs.symbols.gmt"
        ,"c2.cp.reactome.v2023.1.Hs.symbols.gmt" 
        ,"c2.cp.pid.v2023.1.Hs.symbols.gmt"
  ),
  desc=c("GO All"
         ,"GO Biological Process"
         ,"GO Cellular Component"
         ,"GO Molecular Function"
         ,"CP Biocarta"
         ,"CP KEGG"
         ,"CP Reactome"
         ,"CP PID"
  ),
  category=c("C5", "C5", "C5", "C5", "C2", "C2", "C2", "C2")
)# end of data.frame

fwrite(gmt_mtd,'/projectnb/tcwlab/MSigDB/gmt_metadata.csv')

#add pathway/gene set size
pathways<-do.call(rbind,lapply(gmt_mtd$name, function(n){
  paths_list<-fgsea::gmtPathways(file.path('/projectnb/tcwlab/MSigDB/',gmt_mtd[name==n]$gmt))
  paths<-do.call(rbind,lapply(names(paths_list),function(p)data.table(term=p,gene=paths_list[[p]],term.size=length(paths_list[[p]]))))
  return(paths[,source:=n])
  
}))

pathways<-do.call(rbind,lapply(gmt_mtd$name, function(n){
  paths_list<-fgsea::gmtPathways(file.path('/projectnb/tcwlab/MSigDB/',gmt_mtd[name==n]$gmt))
  paths<-do.call(rbind,lapply(names(paths_list),function(p)data.table(term=p,term.size=length(paths_list[[p]]))))
  return(paths[,source:=n])
  
}))

fwrite(pathways,'/projectnb/tcwlab/MSigDB/all_msigdb_pathways.csv.gz')

fwrite(unique(pathways[,.(source,term,term.size)]),'/projectnb/tcwlab/MSigDB/all_msigdb_pathways_metadata.csv.gz')


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
fwrite(genesets_mtd,'/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

fwrite(unique(genesets_mtd[,.(category,subcat,pathway,pathway.size)]),'/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')

