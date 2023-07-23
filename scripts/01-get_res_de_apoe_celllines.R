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

