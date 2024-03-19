source<-function(file,chdir=TRUE)base::source(file,chdir = chdir)

source('r_utils.R')

#UTILITILES####

BamDownSample<-function(bams,out_dir=NULL,prop=0.01){
  bam_dir=unique(dirname(bams))
  if(is.null(out_dir)){
    out_dir=bam_dir
    if(length(out_dir)>1){
      stop('specify outputs directory')
    }
  }
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  outs<-fp(out_dir,str_replace(basename(bams),'.bam$','_head.bam'))
  
  cmds<-paste('samtools view -s',prop, bams,
              '-b -o',outs)
  for(cmd in cmds){
    message('running ',cmd)
    system(cmd)
  }
  message(paste(basename(outs),collapse = ', '),' created in',out_dir)
  
}


#MAIN FUNCTIONS####
#==:> Calculate Splicing Deficiency Score (Sd), based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9833817/ 
#and https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.1310249110&file=pnas.201310249SI.pdf 
# Sd = (# of intronic reads/total intronic length)
# ÷ (# of exonic reads/total exonic length)
#transform to BED regions the bam
#usage: bedtools bamtobed -i file.bam > file.bed'
ConvertBamToBed<-function(bams,out_dir=NULL,job_file=NULL,
                          nThreads=NULL,parallelize=FALSE){
  
  bam_dir=unique(dirname(bams))
  if(is.null(out_dir)){
    out_dir=bam_dir
    if(length(out_dir)>1){
      stop('specify outputs directory')
    }
  }
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  

  cmds<-lapply(bams, function(bam_file){
    paste('bedtools bamtobed -i',bam_file,'|gzip -c >',fp(out_dir,str_replace(basename(bam_file),'.bam$','.bed.gz')))
  })
  
  if(is.null(job_file)){
    
    #job_file<-fp('scripts',ps('bamtobed_',tools::file_path_sans_ext(bams[1],'_andCo.qsub')))
    
    for(cmd in cmds){
      message('running ',cmd)
      system(cmd)
      
    }
    
  }else{
  CreateJobFile(cmds,file = job_file,
                loadBashrc = T,modules = c('bedtools'),nThreads = nThreads,parallelize =parallelize )
  
  jobid<-RunQsub(job_file,job_name = 'bamtobed')
  WaitQsub(job_file,jobid =jobid )
  }
  
  bed_files<-fp(out_dir,str_replace(basename(bams),'.bam$','.bed.gz'))
  
  return(bed_files)
}
#Example:
# astro_files<-list.files('examples_data/',pattern = 'head.bam$',full.names = T)
# bed_files<-ConvertBamToBed(bams = astro_files,out_dir = 'outputs/03-splicing_def_test')
# bed_files<-ConvertBamToBed(bams = astro_files,out_dir = 'outputs/03-splicing_def_test',job_file = 'scripts/03-bamtobed_test.qsub')
# 
# bed_files<-ConvertBamToBed(bams = astro_files[1:3],out_dir = 'outputs/03-splicing_def_test',job_file = 'scripts/03-bamtobed_test.qsub',parallelize = TRUE)

#count number of regions falling in Introns vs Exons canonical rediongs
#canonical regions extract from ucsc
#usage: bedtools intersect -c -wa -a regions.bed -b sample.bed

CountBEDOverlap<-function(bed_files,genomic_regions_file,
                          out_dir=NULL,job_file=NULL,
                          nThreads=NULL,parallelize=F){
  bed_dir=unique(dirname(bed_files))
  if(is.null(out_dir)){
    out_dir=bed_dir
    if(length(out_dir)>1){
      stop('specify outputs directory')
    }
  }
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  
  cmds<-lapply(bed_files, function(b){
    paste('bedtools intersect -c -wa -a',genomic_regions_file,
          '-b',b,'|gzip -c >',
          fp(out_dir,ps(tools::file_path_sans_ext(basename(b)),'.',tools::file_path_sans_ext(basename(genomic_regions_file)),
                        '.overlap.count.bed.gz')))
  })
  
  if(is.null(job_file)){
    
    #job_file<-fp('scripts',ps('counts_overlap_',tools::file_path_sans_ext(bed_files[1],'_andCo_with_',tools::file_path_sans_ext(basename(genomic_regions_file)),'.qsub')))

     for(cmd in cmds){
       message('running ',cmd)
       
       system(cmd)
       
       
     }
    
    }else{
      
    CreateJobFile(cmds,file = job_file,
                  loadBashrc = T,modules = c('bedtools'),nThreads = nThreads,parallelize =parallelize )
    
    jobid<-RunQsub(job_file,job_name = 'countbedoverlap')
    
    WaitQsub(job_file,jobid =jobid )
    }
  
  count_files<- fp(out_dir,ps(tools::file_path_sans_ext(basename(bed_files)),'.',tools::file_path_sans_ext(basename(genomic_regions_file)),
                              '.overlap.count.bed.gz'))
  
  return(count_files)
  


}

#Example:
# introns_counts<-CountBEDOverlap(bed_files =bed_files,
#                                 genomic_regions_file = 'examples_data/some_introns.bed.gz')
# 
# exons_counts<-CountBEDOverlap(bed_files =bed_files,
#                                 genomic_regions_file = 'examples_data/some_exons.bed.gz')
# 

# #overlapping introns and exons
SplicingDeficiency<-function(introns_count,exons_count,out_dir){
  
  # Sd = (# of intronic reads/total intronic length)
  # ÷ (# of exonic reads/total exonic length)
  
 
  counts<-rbind(fread(introns_count,select = c(1:5),col.names = c('chr','start','end','ID','count'))[,type:='intron'],
                fread(exons_count,select = c(1:5),col.names = c('chr','start','end','ID','count'))[,type:='exon'])
  counts[,length:=end-start]
  counts[,count.norm:=count/length]
  counts[,transcript_id:=str_extract(ID,'ENST[0-9]+')]
  counts[,splic.def:=sum(count.norm[type=='intron'])/sum(count.norm[type=='exon']),by=c('transcript_id')]
  
  splic_def<-unique(counts[,.(splic.def,transcript_id)])
  return(splic_def)
  
}

#Example:
# res1<-SplicingDeficiency(introns_counts[1],
#                          exons_count =exons_counts[1])
# 
# 
# res_all<-rbindlist(lapply(1:length(bed_files), function(i){
#   sample<-str_remove(basename(bed_files[i]),'Aligned.sortedByCoord.out_head.bed.gz')
#   
#   
#   res1<-SplicingDeficiency(introns_counts[i],
#                            exons_count =exons_counts[i])[,sample_id:=sample]
#   return(res1)
#   
#   
# }))
# 
# res_all
