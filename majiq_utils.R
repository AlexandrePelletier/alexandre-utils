#MajiqConfig####
#input ; metadata (data.table format) containing bam file paths, and covariate to group by the samples
#output : create <config_file> for `majiq build`
MajiqConfig<-function(mtd_bams,config_file,bam_file_col='bam',
                      group.by='sample_id',genome='hg38',strandness='reverse'){
  bams=mtd_bams[[bam_file_col]]
  groups=mtd_bams[[group.by]]
  bam_dir=unique(dirname(bams))
  if(length(bam_dir)>1){
    stop('bam files should be in the same directory')
  }
  
  bam_ids=tools::file_path_sans_ext(basename(bams))
  
  #create the config file
  cat('[info]\n',file = config_file)
  message('config file created at ',config_file)
  
  #add the general config parameters
  bamdirs_opt<-paste0('bamdirs=',bam_dir)
  genome_opt<-paste0('genome=',genome)
  strandness_opt<-paste0('strandness=',strandness)
  cat( c(bamdirs_opt,genome_opt,strandness_opt),file = config_file,append = T,sep = '\n')
  
  #
  
  #add experiments parameters
  cat('[experiments]\n',file = config_file,append = T)
  bams_groups<-split(bam_ids,groups)
  for(group in names(bams_groups)){
    bams_group=bams_groups[[group]]
    
    cat( paste0(group,'=',paste(bams_group,collapse = ','),'\n'),file = config_file,append = T)
    
    
  }
  
}
'/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid/RawData/organoid/12_S20_L001_R1_001.fastq.gz'

#FormatVoilaHet####
#transform to long format the voila tsv results of majiq heterogen differential splicing
#input: the voila datatable, name of the case/test group and control/reference group
#output: datatable of all splicing event differential results
#note: will also call an `diff_spliced` event if dtPSI>10% and all stats test pvalue<0.05 & 
#`diff_spliced_relaxed` event if dtPSI>5% and at least 3 stats test pvalue<0.05 
FormatVoilaHet<-function(res_voila,case_name,control_name){
  
   
  #first, reformat the PSI columns (containing ';') and take care of harmonizing the dtpsi/HET specific column: 

    res_long<-FormatVoilaDtPsi(res_voila,case_name=case_name,control_name=control_name)

    #add Heterogen specific information    
    res_long[,dt_median_psi:=median_psi_case-median_psi_control]
    res_long[,diff_spliced:=abs(dt_median_psi)>0.1&WILCOXON<0.05&TTEST<0.05&TNOM<0.05]
    res_long[,diff_spliced_relaxed:=abs(dt_median_psi)>0.05&sum(c(WILCOXON<0.05,TTEST<0.05,TNOM<0.05))>1,by='lsv_id']
    
  
  
  return(res_long)
}

FormatVoilaDtPsi<-function(res_voila,case_name,control_name){
  
  #first, reformat the PSI columns (containing ';')
  res_long<-FormatVoilaPSI(res_voila)
  
  #then, take care of harmonizing the dtpsi/HET specific column: 
  cols_control=colnames(res_voila)[str_detect(colnames(res_voila),control_name)]
  cols_case=colnames(res_voila)[str_detect(colnames(res_voila),case_name)&!str_detect(colnames(res_voila),control_name)]
  #cols_groups<-c(cols_control,cols_case)
  
  
  #generalized case and control columns
  for(col in cols_control){
    setnames(res_long,col,paste0(str_remove(col,paste0(control_name,'_')),'_control'))
    
  }
  for(col in cols_case){
    setnames(res_long,col,paste0(str_remove(col,paste0(case_name,'_')),'_case'))
    
  }
  res_long[,case:=case_name]
  res_long[,control:=control_name]
  if('mean_dpsi_per_lsv_junction'%in%colnames(res_voila)){
    res_long[,mean_dpsi_per_lsv_junction:=mean_psi_case-mean_psi_control]
  }
  return(res_long)
}

#FormatVoilaTSV
#transform to long format the voila tsv results of majiq psi/dtpsi/het differential splicing
#input: the voila datatable,
#output: datatable of all splicing event psi quantification
FormatVoilaPSI<-function(res_tsv){
  # cols_stats<-c('TNOM','TTEST','WILCOXON')
  # cols_stats<-c(cols_stats,paste0(cols_stats,'_quantile'))
  cols_tounfold<-setdiff(colnames(res_tsv)[sapply(res_tsv,function(x)any(str_detect(x,';')))],'exons_coords')
  cols_ok<-setdiff(colnames(res_tsv),cols_tounfold)
  
  # cols_tounfold<-c(cols_groups,cols_stats,'tnom_score','changing','nonchanging','de_novo_junctions','junctions_coords')
  
  cols_list<-lapply(cols_tounfold, function(col){
    message(col)
    res_tsv[,strsplit(unlist(.SD),';'),by='lsv_id',.SDcols=col]
    
  })
  
  res_long<-Reduce(cbind,cols_list)
  res_long<-res_long[,.SD,.SDcols=!duplicated(colnames(res_long))]
  #res_long<-merge(res_long,unique(res_tsv[,.(gene_name,gene_id,lsv_id,lsv_type,num_junctions,num_exons,seqid,strand,exons_coords,ir_coords,ucsc_lsv_link)]))
  #add the remaining cols
  res_long<-merge(res_long,unique(res_tsv[,.SD,.SDcols=cols_ok]))
  
  #transform numerical columns 
  cols_num<-c(colnames(res_long)[str_detect(colnames(res_long),'mean|median|quantile|percentile|TNOM|TTEST|WILCOXON|psi|probability')],'de_novo_junctions')
  res_long[,(cols_num):=lapply(.SD,as.numeric),.SDcols=cols_num]
  return(res_long)
}


#AnnotVoila####
#annotate with modulizer results
#Inputs: data.table of long format of voila results (pass through `FormatVoila` functions)
#outputs: annotated data.table, can have more row than the input data,table, because one junction can be annotated with 2 different splicing_type depending of the context
AnnotVoila<-function(voila_long,modulizer_dir,annos=c('cassette','intron','alt3','alt5','first','last','orphan','exclusive'),
                     cols_to_keep=c('lsv_id','complex','denovo','spliced_with','spliced_with_coord','junction_name','junction_coord','event_size','event_type')){
  anno_files<-list.files(modulizer_dir,pattern = paste(annos,collapse = '|'),full.names = T)
  anno<-rbindlist(lapply(anno_files,
                         function(f)fread(f)[,event_type:=tools::file_path_sans_ext(basename(f))]),fill = T)
  anno<-unique(anno)[lsv_id!='']
  
  if(!is.null(cols_to_keep)){
    anno<-anno[,.SD,.SDcols=cols_to_keep]
    
  }
  voila_long_anno<-merge(voila_long,anno,all.x = T,allow.cartesian = T)
  voila_long_anno[is.na(event_type),event_type:='other']
  voila_long_anno[,junction_type:=ifelse(de_novo_junctions==1,'denovo','known')]
  voila_long_anno[,junction_type:=factor(junction_type,levels=c('known','denovo'))]
  return(voila_long_anno)
}
