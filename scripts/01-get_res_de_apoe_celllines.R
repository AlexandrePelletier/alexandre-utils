out<-'outputs/01-res_de_apoe_cell_lines'
dir.create(out)
source('r_utils.R')
iso_path<-'/projectnb/tcwlab/Project/IsogenicAPOE/Analysis/DESeq2/'
celltypes<-c('Astro','Micro','Neuro_composition')
res_de<-Reduce(rbind,lapply(celltypes,function(ct){
  indivs<-list.files(fp(iso_path,ct))
  
  res_ct<-Reduce(rbind,lapply(indivs,function(ind){
    res_dir<-fp(iso_path,ct,ind)
    res<-fread(fp(res_dir,list.files(res_dir,pattern = 'NoGeneIDDuplicates')))
    return(res[,indiv:=ind][,comparison:=ifelse(str_detect(ind,'KO'),'KOvs33','44vs33')])
  }))
  return(res_ct[,cell_type:=ct])
}))
res_de[comparison=='KOvs33']
res_de[cell_type=='Micro'][indiv=='Indiv1']
fwrite(res_de,fp(out,'res_deseq2_APOE_isogenics.csv'))
res_de<-fread(fp(out,'res_deseq2_APOE_isogenics.csv'))
table(res_de[padj<0.05&abs(log2FoldChange)>0.5][,.(indiv,cell_type)])

pop_path<-'/projectnb/tcwlab/Project/PopulationAPOE/JuliaNote/Analysis/Deseq2_GRSnoAPOE/'
celltypes<-c('Astro','Micro','Neuro_composition','Neuro')
res_de_pop<-Reduce(rbind,lapply(celltypes,function(ct){
  indivs<-list.files(fp(iso_path,ct))
  
    res_dir<-fp(pop_path,ct)
    res_ct<-fread(fp(res_dir,list.files(res_dir,pattern = 'NoGeneIDDuplicates')))

  return(res_ct[,cell_type:=ct])
}))
res_de_pop[,comparison:='44vs33']

fwrite(res_de_pop,fp(out,'res_deseq2_APOE_population_corrGRSNoAPOE.csv'))




pop_path<-'/projectnb/tcwlab/Project/PopulationAPOE/JuliaNote/Analysis/Deseq2_GRSnoAPOE/'
celltypes<-c('Astro','Micro','Neuro_composition','Neuro')
res_de_pop<-Reduce(rbind,lapply(celltypes,function(ct){
  indivs<-list.files(fp(iso_path,ct))
  
  res_dir<-fp(pop_path,ct)
  res_ct<-fread(fp(res_dir,list.files(res_dir,pattern = 'NoGeneIDDuplicates')))
  
  return(res_ct[,cell_type:=ct])
}))
res_de_pop[,comparison:='44vs33']

fwrite(res_de_pop,fp(out,'res_deseq2_APOE_population_corrGRSNoAPOE.csv'))


#mouseApoe
dirs_path<-'/projectnb/tcwlab/Project/mouseApoe/Analysis/DESeq2/'
deseq2_dirs<-list.dirs(dirs_path,full.names = T,recursive = F)
celltypes<-c('Astro','Micro','Neuro_composition','Neuro')
res_de_mouse<-Reduce(rbind,lapply(deseq2_dirs,function(dir){
  cond=basename(dir)
  res_cond<-fread(fp(dir,list.files(dir,pattern = 'NoGeneIDDuplicates')))
  
  return(res_cond[,condition:=cond])
}))

res_de_mouse<-res_de_pop

fwrite(res_de_mouse,fp(out,'res_deseq2_APOE_Mouse.csv'))


#scRNAseq
resde<-fread('../projects/juao_autophagy/outputs/13-isogenics_analysis/res_sc_iso_44vs33_tcw1to3.csv.gz')
table(resde$cell_type)

fwrite(resde,fp(out,'res_deseq2_singlecell_MCC_isogenics_tcw1and2merged.csv'))

#check some genes
res_de<-fread(fp(out,'res_deseq2_APOE_isogenics.csv'))
res_de[padj<0.1]

res_de[gene_name%in%c('SERPINA3','LCN2','NGAL')]
ENSG00000148346
res_de[str_extract(V1,'ENSG[0-9]+')=='ENSG00000148346']




