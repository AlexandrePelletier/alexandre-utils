
PrepDESeq2Design<-function(mtd,covs,scale_num=TRUE){
  
  #remove if all sample have the same value
  to_keep<-sapply(covs,function(cov)length(unique(mtd[[cov]]))!=1)
  if(sum(!to_keep)>1)message('remove ',paste(covs[!to_keep],collapse = ', '),'because all samples have the same value')
  covs_filtered<-covs[to_keep]
  
  #if a factor, need to have at least 2 sample in each category
  to_keep<-sapply(covs_filtered,function(cov)ifelse(!is.numeric(mtd[[cov]]),all(table(mtd[[cov]])!=1),TRUE))
  if(sum(!to_keep)>1)message('remove ',paste(covs_filtered[!to_keep],collapse = ', '),' factor(s) because need to have at least 2 samples in each category')
  
  covs_filtered<-covs_filtered[to_keep]
  
  #if numerical, scales the data
  
  if(scale_num){
    numerical_factors<-covs_filtered[mtd[,sapply(.SD,is.numeric),.SDcols=covs_filtered]]
    
    if(length(numerical_factors)>0){
      numerical_factors_scaled<-paste0(numerical_factors,'_scaled')
      message('scaling ',paste(numerical_factors,collapse=', '),' in new column ',paste(numerical_factors_scaled,collapse=', '))
      mtd[,(numerical_factors_scaled):=lapply(.SD,scale),.SDcols=numerical_factors]
      covs_filtered<-sapply(covs_filtered, function(x)ifelse(x%in%numerical_factors,paste0(x,'_scaled'),x))
    }
    
  }
 
  formul_clean<-formula(paste0("~", paste(covs_filtered,collapse ='+')))
  return(formul_clean)

  
}

ScaleDESeq2Covs<-function(mtd,covs){
  mtd_scaled<-copy(mtd)
  numerical_factors<-covs[mtd_scaled[,sapply(.SD,is.numeric),.SDcols=covs]]
  numerical_factors_scaled<-paste0(numerical_factors,'_scaled')
  message('scaling ',paste(numerical_factors,collapse=', '),' in new column ',paste(numerical_factors_scaled,collapse=', '))
  mtd_scaled[,(numerical_factors_scaled):=lapply(.SD,scale),.SDcols=numerical_factors]
  return(mtd_scaled)
}


#RunFgseaMsigdb
#inputs: res_de: datatable of differential expression results, with gene name (in gene symbol format) and statistical test
#       msigdb_path : path to MSigdb data, should a table with at least column 'pathway', 'gene', and 'category' of the pathways
#       score : from which column of the res_de to collect the score/statistic to use to rank genes for fgsea?
#outputs : res of fgsea in data table format
#Notes : fgsea will be run by 'category' of the Msigdb (Canonical Pathways 'CP' and Gene ontoloy 'GO'  term by default) if  can specify 
#duplicate_choice: how deal with duplicated gene_names / id/.  'random'ly pick one (default), or keep only the 'top' one based on the stat, or do 'nothing'. 
RunFgseaMsigdb<-function(res_de,score='stat',rankbased=F,
                         msigdb_path='/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz',
                         genes_cols=c('gene','gene_name','gene_id'),
                         duplicate_choice='random',
                         group.by=NULL,minSize = 10,maxSize = 2000,
                         gseaParam = 1,scoreType='std',eps=1e-50,
                         nPermSimple = 10000,...){
  require(data.table)
  require(fgsea)
  if(!'data.table'%in%class(res_de)){
    if(!any(genes_cols%in%colnames(res_de))){
      res_de$gene<-rownames(res_de)
    }
    res_de<-data.table(res_de)
  }else{
    res_de<-copy(res_de)
  }
  
  gene_col=which(colnames(res_de)%in%genes_cols)[1]
  res_de$gene<-res_de[[gene_col]]
  
  if(!is.null(group.by)){
   
    res_de_list<-split(res_de,by=group.by)
    
    res_fgsea<-rbindlist(lapply(names(res_de_list),function(g){
      message('testing msigdb pathway enrichment in ',g)
      res_gsea<-RunFgseaMsigdb(res_de_list[[g]],score=score,
                               msigdb_path=msigdb_path,
                               genes_cols=genes_cols)
      return(res_gsea[,query:=g])
      
    }))
    
  }else{
    msigdb<-fread(msigdb_path)
    
    to_rm<-is.na(res_de[[score]])
    if(sum(to_rm)>0){
      warning('removing ',sum(to_rm),' genes containing missing ',score,' value.')
      res_de<-res_de[!to_rm]
      
    }
    dups<-duplicated(res_de[order(-res_de[[score]])]$gene_name)
    
    if(sum(dups)>0){
      if(sum(dups)>0.75*length(unique(res_de$gene_name))){
        stop('>75% of duplicated gene names, probably several non separated conditions.')
      }
      warning(sum(dups),' duplicated genes')
      if(duplicate_choice=='top'){
        warning(' removing them by picking the top stat per gene')
        res_de<-res_de[order(-res_de[[score]])][!(dups)]
        
        
      }else if(duplicate_choice=='random'){
        warning(' removing them by randomly picking one unique gene')
        res_de<-unique(res_de[sample(1:nrow(res_de))],by=gene_col)
      
      }
    }
    
    if(rankbased){
      stats<-setNames(sign(res_de[[score]])*rank(abs(res_de[[score]])),res_de$gene)
    }else{
      stats<-setNames(res_de[[score]],res_de$gene)
      print(head(sort(stats,decreasing = T)))
      
    }
    
    res_fgsea<-rbindlist(lapply(unique(msigdb$category), function(cat){
      
      msigdbf<-msigdb[category==cat]
      pathways=split(msigdbf$gene,msigdbf$pathway)
      res<-fgsea(pathways,stats,minSize = minSize,maxSize = maxSize,
                 gseaParam = gseaParam,scoreType=scoreType,eps=eps,
                 nPermSimple = nPermSimple,...)
      
      return(res[,category:=cat])
    }))
    
    #annot 
    res_fgsea<-merge(res_fgsea,unique(msigdb[,-'gene']))[order(padj)]
    
  }
  
  
  return(res_fgsea)
}
