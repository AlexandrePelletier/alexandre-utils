
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
  vecx<-strsplit(x,sep)[[1]]
  vecy<-strsplit(y,sep)[[1]]
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
