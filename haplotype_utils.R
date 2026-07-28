

#getLocalHaplo: from topmed phased vcf.gz file, extract a region haplotype,returning a table the haplotype ID and SNPs for the 2 chromosomes of each sample
#default QC: AF >2.5% & imputed R2>0.8
# 
getLocalHaplo<-function(vcf,region,samples_mtd,indiv_col='IID',min_AF=0.025,min_imput_R2=0.8){
  
  geno<-fread(paste('tabix -h',vcf,region,
                    '|',"grep -v '##'"))
  #QC
  geno[,maf:=as.numeric(str_extract(INFO,'[0-9.e-]+'))]
  geno[,R2:=as.numeric(str_remove(str_extract(INFO,'R2=[0-9.]+'),'R2='))]
  geno[,genotyped:=str_detect(INFO,'ER2')]
  message('Filtering SNPs: maf: ',round(min_AF*100,digits = 1),'%',', R2 imputation: ',min_imput_R2)
  
  genof<-geno[maf>min_AF&R2>min_imput_R2]
  message(nrow(genof),' SNPs passing QC')
  
  genolong<-merge(genof[,1:8],melt(genof[,.SD,.SDcols = c(samples_mtd[[indiv_col]],'ID')],id.vars = 'ID',variable.name = indiv_col,value.name = 'genotype_info'))
  genolong[,dose:=str_extract(genotype_info,'[01]\\|[01]')|>str_count('1')]
  genolong<-merge(genolong,samples_mtd,by=indiv_col)
  
  genolong[,genotype:=strsplit(genotype_info,'\\:')[[1]][1],by=c(indiv_col,'POS','ID')]
  
  #then separate allele 1 and 2
  genolong[,A1:=as.numeric(strsplit(genotype,'\\|')[[1]][1]),by=c(indiv_col,'POS','ID')]
  genolong[,A2:=as.numeric(strsplit(genotype,'\\|')[[1]][2]),by=c(indiv_col,'POS','ID')]  
  
  homologs<-rbind(setnames(genolong[order(POS,ID),ifelse(A1==0,REF,ALT),by=c(indiv_col,'POS','ID')],'V1','SNP')[,haplot:='H1'],
                  setnames(genolong[order(POS,ID),ifelse(A2==0,REF,ALT),by=c(indiv_col,'POS','ID')],'V1','SNP')[,haplot:='H2'])
  homologs[,homolog_id:=paste(.SD[[1]],haplot,sep = '_'),.SDcols=indiv_col]
  
  # get an unique haplotype ID 
  haplots<-rbind(genolong[order(POS,ID),list(haplo=paste(ifelse(A1==0,REF,ALT),collapse = '_')),by=indiv_col][,haplot:='H1'],
                 genolong[order(POS,ID),list(haplo=paste(ifelse(A2==0,REF,ALT),collapse = '_')),by=indiv_col][,haplot:='H2'])
  
  haplo_ids=data.table(haplo=unique(haplots$haplo),
                       haplo_id=paste0('haplo',1:length(unique(haplots$haplo))))
  haplo_dt<-merge(haplots,haplo_ids)
  haplo_dt
  
  homologs<-merge(homologs,haplo_dt[,-'haplo'])
  
  #add sample anno
  
  homologs<-merge(homologs,samples_mtd,by=indiv_col)
  
  #add SNP anno
  homologs<-merge(homologs,genof[,1:8],by=c('ID','POS'))
  homologs[,maf:=as.numeric(str_extract(INFO,'[0-9.e-]+'))]
  homologs[,R2:=as.numeric(str_remove(str_extract(INFO,'R2=[0-9.]+'),'R2='))]
  homologs[,genotyped:=str_detect(INFO,'ER2')]
  
  return(homologs)
  
}

#getHaplo
#get Haplotype from phased vcf table
#Arguments:
# vcf: the phased vcf to extract haplotype for all samples
# haplo_ref: if a haplotype have already been produced for reference, will annotate with the haplo_id of this reference if possible. 

#return long table where each row is a SNP from an haplotype of a sample annotated from an unique haplo id
GetHaplo<-function(vcf,haplo_ref=NULL){
  nt_codes<-data.table(SNP=c("A","T","C","G"),
                       snp_code=1:4)
  samples<-setdiff(colnames(vcf),c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','maf','R2','genotyped'))
  
  
  haplot_dt<-melt(vcf,measure.vars = samples,variable.name = 'IID',value.name ='genotype_info')
  haplot_dt 
  #from the full genotype and QC info (eg : 0|0:1,0,0:0:0,0), extract the genotype only (first element)
  haplot_dt[,genotype:=strsplit(genotype_info,'\\:')[[1]][1],by=c('IID','POS','ID')]
  
  #then separate haplo 1 and 2
  haplot_dt[,H1:=as.numeric(strsplit(genotype,'\\|')[[1]][1]),by=c('IID','POS','ID')]
  haplot_dt[,H2:=as.numeric(strsplit(genotype,'\\|')[[1]][2]),by=c('IID','POS','ID')]

  #long the format again to have one snp-sample-haplo by line. Allele being the nucleotide that the sample have at each position
  snps_haplos<-rbind(haplot_dt[order(POS,ID),.(SNP=ifelse(H1==0,REF,ALT)),by=c('IID','POS','ID')][,haplo:='H1'],
                   haplot_dt[order(POS,ID),.(SNP=ifelse(H2==0,REF,ALT)),by=c('IID','POS','ID')][,haplo:='H2'])
  
  # get an unique haplotype ID 

    haplo_ext<-rbind(haplot_dt[order(POS,ID),list(haplo_seq=paste(ifelse(H1==0,REF,ALT),collapse = '_')),by='IID'][,haplo:='H1'],
                     haplot_dt[order(POS,ID),list(haplo_seq=paste(ifelse(H2==0,REF,ALT),collapse = '_')),by='IID'][,haplo:='H2'])
    message(length(unique(haplo_ext$haplo_seq)),' unique haplotype found')
    
    haplo_ids=data.table(haplo_seq=unique(haplo_ext$haplo_seq),
                         haplo_id=paste0('haplo',1:length(unique(haplo_ext$haplo_seq))))
    haplo_dt<-merge(haplo_ext,haplo_ids)
    snps_haplos<-merge(snps_haplos,haplo_dt)
    
    if(!is.null(haplo_ref)){
      #filter for SNP present in ref and get the haplo sequence for the query sample
      snps_haplosf<-merge(snps_haplos[,.(IID,POS,SNP,haplo)],
                             unique(haplo_ref[,.(POS,SNP)]),by=c('POS','SNP'))[order(POS,SNP)]

      snps_haplosf[,haplo_seq:=paste(SNP,collapse = '_'),by=c('IID','haplo')]
      #find if this haplo seq is found in the ref
      haplo_ref[,haplo_id.ref:=haplo_id]
      
     haplosf<- merge(unique(snps_haplosf[,.(IID,haplo,haplo_seq)]),
                           unique(haplo_ref[,.(haplo_seq,haplo_id.ref)]),by=c('haplo_seq'))[,-'haplo_seq']
      
      message('was able to find match for ',length(unique(haplosf$haplo_id.ref)),
              '/',length(unique(haplo_ref$haplo_id)),' ref haplotypes')
      
      if(length(unique(haplosf$haplo_id.ref))>0){
        #add the haplo_ref id to the query 
        snps_haplos<-merge(snps_haplos,haplosf,by=c('IID','haplo'),all.x=TRUE)
      }
    }
    
  #create sample-haplotype ID
  snps_haplos[,sample_haplo:=paste(IID,haplo,sep='_')]

  return(snps_haplos)
  
}


#DNA sequence analysis#### 
#hamming_distance
#cmpute hamming distance (number of mismatch) for 2 string sequence of same size 
hamming_distance<-function(x,y,sep=''){
  if(length(x)==1){
    vecx<-strsplit(x,sep)[[1]]
    vecy<-strsplit(y,sep)[[1]]
  }else{
    vecx=x
    vecy=y
  }

  sum(vecx!=vecy)
}
#if consider indels, can use getDist of Alakazam
#install.packages('alakazam')
#need format seq like: ATC.G with '.' being gap
# alakazam::seqDist(seq1, seq2, dist_mat = getDNAMatrix(gap=-1))#add +1 per base insertion

# AlignSeqs<-function(seq_with_underline1,seq_with_underline2){
#   vec<-strsplit(seq_with_underline,'_')[[1]]
#   return(seq_with_dotgap)
# }



heatmapHaplo<-function(x,variant_col='ID',annotation_indiv=NULL,annotation_snp=NULL,display_SNP=NULL,cluster_rows=TRUE){
  
  nt_codes<-data.table(SNP=c("A","T","C","G"),
                       snp_code=1:4)
  
  x<-merge(x,nt_codes,by='SNP')
  
  mtd_haplo<-unique(x,by='homolog_id')
  toconv<-colnames(mtd_haplo)[sapply(mtd_haplo, is.logical)]
  
  if(length(toconv)>0){
    setnames(mtd_haplo,toconv,paste0(toconv,'1'))
    mtd_haplo[,(toconv):=lapply(.SD,as.numeric),.SDcols=paste0(toconv,'1')]    
  }
  
  mtd_snps<-unique(x,by='ID')
  toconv<-colnames(mtd_snps)[sapply(mtd_snps, is.logical)]
  if(length(toconv)>0){
    setnames(mtd_snps,toconv,paste0(toconv,'1'))
    mtd_snps[,(toconv):=lapply(.SD,as.numeric),.SDcols=paste0(toconv,'1')]
    
  }
  if(cluster_rows){
    hamm_dist<-sapply(unique(x$homolog_id),function(h1){
      sapply(unique(x$homolog_id), function(h2){
        hamming_distance(x[h1,on='homolog_id'][order(POS)]$SNP,x[h2,on='homolog_id'][order(POS)]$SNP)
        
      })
    })
    hamm_dist<-as.dist(hamm_dist)
    
  }else{
    hamm_dist='euclidean'
  }

  
  haplo_map<-dcast(x,homolog_id~ID,value.var = 'snp_code')
  haplo_snp<-dcast(x,homolog_id~ID,value.var = 'SNP')
  
  snps<-colnames(haplo_map)
  if(is.null(display_SNP)){
    display_SNP=length(snps)<60
  }
  if(display_SNP){
    snpmat=data.frame(haplo_snp,row.names = 'homolog_id')[mtd_haplo$homolog_id,mtd_snps[order(POS)]$ID]
  }else{
    snpmat=FALSE 
  }
  
  if(!is.null(annotation_indiv)){
    
    anno_homol=data.frame(mtd_haplo,row.names = 'homolog_id')[,annotation_indiv,drop=FALSE]
  }else{
    anno_homol=NA
  }
  
  if(!is.null(annotation_snp)){
    anno_snp=data.frame(mtd_snps,row.names = 'ID')[,annotation_snp,drop=FALSE]
  }else{
    anno_snp=NA
  }
  
  pheatmap(data.frame(haplo_map,row.names = 'homolog_id')[mtd_haplo$homolog_id,mtd_snps[order(POS)]$ID],
           color = c('green3','cyan2','orange3','coral1'),
           annotation_row =anno_homol,
           annotation_col =anno_snp,
           cluster_cols  = F,border_color = 'grey',
           cluster_rows = cluster_rows,
           clustering_distance_rows=hamm_dist,
           legend = TRUE, legend_breaks = c(1.5,2.25,3,3.75),legend_labels = c('A','T','C','G'),
           show_colnames = T,
           display_numbers =snpmat
  )
  
  
}
