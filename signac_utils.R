source<-function(file,chdir=TRUE)base::source(file,chdir = chdir)

source('r_utils.R')

####Signac####


CheckMotif<-function(object,peaks,motif.name=NULL,assay = NULL,return.peaks=FALSE){
  require("Signac")
  if(is.null(assay))assay<-DefaultAssay(object)
  motif<-GetMotifIDs(object,motif.name,assay=assay)
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data")
  
  motifs_peaks_tf <- motif.all[peaks,motif , drop = FALSE]
  
  if(return.peaks){
    motifs_peaks_tf<-rownames(motifs_peaks_tf)[as.vector(motifs_peaks_tf==1)]
    return(motifs_peaks_tf)
  }else{
    motifs_peaks_tf_vec<-as.vector(motifs_peaks_tf==1)
    names(motifs_peaks_tf_vec)<-rownames(motifs_peaks_tf)
    return(motifs_peaks_tf_vec)
  }
  
  
}
GetMotifIDs<-function(object,motif.names,assay=NULL,return_dt=FALSE){
  if(is.null(assay))assay<-DefaultAssay(object)
  idx<-match(motif.names,object@assays[[assay]]@motifs@motif.names)
  if(return_dt){
    return(
      data.table(motif.name=motif.names,
                 motif.id=names(object@assays[[assay]]@motifs@motif.names[idx]))
    )
  }else{
    return(names(object@assays[[assay]]@motifs@motif.names[idx]))
  }
  
}

GetMotifs<-function(object,peaks,motifs=NULL,assay = NULL){
  #return tfmotif peak data.table
  require("Signac")
  require("data.table")
  
  if(is.null(assay))assay=DefaultAssay(object)
  
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  
  if(is.null(motifs)) motifs<-colnames(motif.all)
  
  motif.filtered<-motif.all[peaks,motifs,drop=F]
  
  motif_dt<-melt(data.table(as.matrix(motif.filtered),keep.rownames = "peak"),
                 variable.name="motif",value.name = "presence")
  
  motif_dt[,presence:=as.logical(presence)]
  motif_dt[(presence)] 
  motif_dt<-motif_dt[(presence)][,-"presence"]
  
  motifsnames<-data.table(motif.name=object[[assay]]@motifs@motif.names,
                          motif=names(object[[assay]]@motifs@motif.names))
  motif_dt<-merge(motif_dt,motifsnames)
  
  ranges_list<-object[[assay]]@motifs@positions[ unique(motif_dt$motif)]

  motifs_pos <- rbindlist(lapply(names(ranges_list),function(motifid){
    
    motif_pos<-rbindlist(lapply(peaks, function(peak){
      start.pos <- start(peak)
      end.pos <- end(peak)
      chromosome <- seqid(peak)
      data.table(as.data.frame(ranges_list[[motifid]]))[seqnames==chromosome&start>start.pos&end<end.pos][,motif:=motifid][,peak:=peak]
      
    }))
    
    
  }))
  
  motif_dt=merge(motif_dt,motifs_pos,by=c('motif','peak'))

  
  
  return(motif_dt)
}








TFMotifPlot<-function(object,region,motif.name,assay=NULL){
  if(is.null(assay))assay<-DefaultAssay(object)
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  ranges<-object@assays[[assay]]@motifs@positions[[GetMotifIDs(object,motif.names = motif.name)]]
  dt <- data.table(as.data.frame(ranges))
  dt_reg<-dt[seqnames==chromosome&start>start.pos&end<end.pos]
  
  p<-ggplot(data = dt_reg) + geom_segment(aes(x = start, y = 0, 
                                              xend = end, yend = 0),col="black", size = 2, data = dt_reg)
  
  p<-p+ theme_classic() + ylab(label = motif.name) + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  return(p)
}

TFsMotifPlot<-function(object,region,motif.names,assay=NULL,size=2,alpha=1,pad=0,expand=FALSE){
  if(is.null(assay))assay<-DefaultAssay(object)
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  
  dt_region <-Reduce(rbind,lapply(motif.names,function(x){
    
    ranges<-object@assays[[assay]]@motifs@positions[GetMotifIDs(object,motif.names = x)]
    
    dt<-data.table(as.data.frame(ranges))[seqnames==chromosome&start>start.pos&end<end.pos][,motif.name:=x]
    return(dt)
  }
  ))
  if(expand){
    p<-ggplot(data = dt_region) + geom_segment(aes(x = start-pad, y = as.numeric(as.factor(motif.name)), 
                                                   xend = end+pad, yend = as.numeric(as.factor(motif.name)),col=motif.name), size = size,alpha=alpha, data = dt_region)
    
  }else{
    p<-ggplot(data = dt_region) + geom_segment(aes(x = start-pad, y = 0, 
                                                   xend = end+pad, yend = 0,col=motif.name), size = size,alpha=alpha, data = dt_region)
    
  }
  
  p<-p+ theme_classic() + ylab(label = "TF motif") + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  return(p)
}



#ExtractCellsFrag: from a tab delimited .gz fragment file, extract fragment of a specific cell
#input: your fragment file, the cells to extract
#output: a <fragment_file>.<name>.extract.tsv.gz file containing only the fragment for this paritcular set of cells,
#if return_file=TRUE,  the output file will be load in R as a data.table, else, only the output file path is return

#output:
ExtractCellsFrag<-function(file,cells,extraction.name='cells',out_dir='.',return_file=FALSE,clean_out_file=FALSE){
  
  out_fragment_file<-file.path(out_dir,paste0(tools::file_path_sans_ext(basename(file),compression = T),'.',extraction.name,'.tsv.gz'))
  
  cmd<-paste('zcat',file,'|grep -Ew',paste0('"',paste(cells,collapse = '|'),'"'),
             '|gzip -c',
             '>',out_fragment_file)
  system(cmd)
  
  if(return_file){
    cells_fragment<-fread(out_fragment_file)
    if(clean_out_file){
      system(paste('rm',out_fragment_file))
    }
    return(cells_fragment)
  }else{
    return(out_fragment_file)
  }
}
#test for astrocytes
# ExtractCellsFrag(frag@path,cells = frag@cells)
# fread('D19-12513.fragments.cells.tsv.gz')


#CountFragmentBedOverlap
#inputs : 1) seurat object contaning fragments objects 
#         2) bed files of the genomics regions to count overlap
#for each genomic region, count how many read overlap
#outputs : bed file <bed_file_name>.<genomics_regions_name>.overlap.count.bed.gz of the number of read falling in each genomic region

#Options:
#         if reduce =TRUE will output only a vector of the total read count overlapping the genomics region for each cell

CountFragmentsBedOverlap<-function(x,genomic_regions_file,name='cells',
                                   reduce=TRUE,
                                   group.by.genomicID=FALSE,
                                   tmp_dir='tmp',
                                   clean_tmp_dir=TRUE,
                                   parallelize=F,
                                   job_file='X-count_fragments_overlap.qsub',
                                   job_name='fragOverlap',
                                   nThreads=NULL,
                                   maxChildJobs=40,
                                   wait_job=TRUE,
                                   wait_for=NULL){
  if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir)
  }
  
  if('Seurat'%in%class(x)){
    if(length(Fragments(x))==0){
      stop('the Seurat object default assay do not contain the fragments objects')
    }
    fragments=Fragments(x)
  }else{
    stop('input file is not a Seurat object')
  }
  
  #create the fragment bed file
  
  fragfiles<-sapply(1:length(fragments), function(i){
    frag=fragments[[i]]
    
    sub_fragfile<-ExtractCellsFrag(file = frag@path,
                     cells =frag@cells,
                     extraction.name = name,
                     out_dir = tmp_dir,return_file = F)
    
    
    return(sub_fragfile)
    
    
  })
  
  #count overlap
  overlap_files<-CountBEDOverlap(bed_files = sub_fragfile,
                  genomic_regions_file = genomic_regions_file,
                  write_bed_file_intervals = TRUE, 
                  out_dir = tmp_dir,
                  parallelize = parallelize,
                  job_file = job_file,
                  job_name = job_name,
                  nThreads = nThreads,
                  maxChildJobs = maxChildJobs,
                  wait_job = wait_job,
                  wait_for = wait_for
  )
 res_counts_list<-lapply(overlap_files, function(file){
    res_counts<-fread(file,select = c(1:4,8,9),col.names = c('chr','start','end','ID','cell_id','count'))
    if(reduce){
      if(group.by.genomicID){
        res_counts<-res_counts[,sum(count),by=c('cell_id','ID')]
        
      }else{
        res_counts<-res_counts[,sum(count),by='cell_id']
        
      }
    }
    return(res_counts)
    })
 
 #get seurat cell id
 res_counts_list<-lapply(1:length(fragments), function(i){
   frag=fragments[[i]]
   res_counts=res_counts_list[[i]]
   
   res_counts[,cell_id_seurat:=names(frag@cells)[match(cell_id,frag@cells)]]
   
   return(res_counts)
   
   })
 
 res_counts_list<-rbindlist(res_counts_list)
 
 #clean tmp dir
 if(clean_tmp_dir){
   system(paste('rm -rf',tmp_dir))
 }
 return(res_counts_list)
  
  
}
