library(data.table)
library(fgsea)

gmt_mtd <- data.table(
  name=c("GO_all", "GO_bp", "GO_cc", "GO_mf",'CP_all', "CP_biocarta", "CP_kegg", "CP_reactome", "CP_pid"),
  
  
  gmt=c("c5.go.v2023.1.Hs.symbols.gmt"
        ,"c5.go.bp.v2023.1.Hs.symbols.gmt" 
        ,"c5.go.cc.v2023.1.Hs.symbols.gmt" 
        ,"c5.go.mf.v2023.1.Hs.symbols.gmt"
        ,'c2.cp.v2023.1.Hs.symbols.gmt'
        ,"c2.cp.biocarta.v2023.1.Hs.symbols.gmt"
        ,"c2.cp.kegg.v2023.1.Hs.symbols.gmt"
        ,"c2.cp.reactome.v2023.1.Hs.symbols.gmt" 
        ,"c2.cp.pid.v2023.1.Hs.symbols.gmt"
  ),
  desc=c("GO All"
         ,"GO Biological Process"
         ,"GO Cellular Component"
         ,"GO Molecular Function"
         ,'CP All'
         ,"CP Biocarta"
         ,"CP KEGG"
         ,"CP Reactome"
         ,"CP PID"
  ),
  category=c("C5", "C5", "C5", "C5", 'C2',"C2", "C2", "C2", "C2")
)# end of data.frame

fwrite(gmt_mtd,'/projectnb/tcwlab/MSigDB/gmt_metadata.csv')

#add pathway/gene set size
pathways<-do.call(rbind,lapply(gmt_mtd$name, function(n){
  paths_list<-fgsea::gmtPathways(file.path('/projectnb/tcwlab/MSigDB/',gmt_mtd[name==n]$gmt))
  paths<-do.call(rbind,lapply(names(paths_list),function(p)data.table(term=p,gene=paths_list[[p]],term.size=length(paths_list[[p]]))))
  return(paths[,source:=n])
  
}))


fwrite(pathways,'/projectnb/tcwlab/MSigDB/all_msigdb_pathways.csv.gz')

fwrite(unique(pathways[,.(source,term,term.size)]),'/projectnb/tcwlab/MSigDB/all_msigdb_pathways_metadata.csv.gz')

