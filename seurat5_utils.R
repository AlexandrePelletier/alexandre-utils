#SeuratV5 utils
require(Seurat)


CalcSplicDef <- function(obj, ...) {
  UseMethod("CalcSplicDef")
}

#CalcSplicDef.Seurat
#calculate the splicing ratio as a new assay and create metadata with splicing_deficiency score.
#this score is calculate per group of cells (eg celltype) which should share similar transcriptome profile
#obj: seurat objecti
#unspliced.assay: name of the unpliced count assay
#spliced.assay: name of the spliced count assay
#allcount.assay: names of the allcount.assay if spliced assay not present
#group.by: per which metadata annotation to calculate the splicing deficit score
#
CalcSplicDef.Seurat<-function(obj,unspliced.assay='unspliced',
                              spliced.assay=NULL,
                              allcount.assay=NULL,
                              group.by='cell_type',
                              min.cells=20,min.genes=50,
                              split.by=NULL,on.disk=NULL,
                              bpcells_dir=NULL,force_bpcells=FALSE,
                              nThreads=1,
                              assay.name='splicingdef',
                              silent=FALSE,
                              ...){
  require(parallel)
  require(BPCells)
  
  tot.cells<-ncol(obj)
  
  with_bp_cells<-tot.cells>200000
  
  if(is.null(on.disk))
    on.disk=FALSE
    
  if(with_bp_cells|on.disk==TRUE){
    on.disk=TRUE
    if(is.null(bpcells_dir))stop('Save the splicing deficit matrix on disk, need the bpcells_dir parameter')
    
    if(!dir.exists(bpcells_dir)){
      dir.create(bpcells_dir)
    }
    
    }
  
  #if split
  if(!is.null(split.by)){
    message('Splitting data per ',split.by)
    
    if(on.disk){
      #work with layers : split assays per layer according to the split.by arg if not already the case
      if(!all(unique(unlist(obj[[split.by]]))%in%str_remove(Layers(obj[[unspliced.assay]],search='counts'),'counts\\.'))){
        if(length(Layers(obj[[unspliced.assay]],search='counts'))>1)
          obj[[unspliced.assay]]<-JoinLayers(obj[[unspliced.assay]])
        
        obj[[unspliced.assay]]<-split(obj[[unspliced.assay]],unlist(obj[[split.by]]))
        
      }
      
      if(!is.null(spliced.assay)){
        
        if(!all(unique(unlist(obj[[split.by]]))%in%str_remove(Layers(obj[[unspliced.assay]],search='counts'),'counts\\.'))){
          if(length(Layers(obj[[spliced.assay]],search='counts'))>1)
            obj[[spliced.assay]]<-JoinLayers(obj[[spliced.assay]])
          obj[[spliced.assay]]<-split(obj[[spliced.assay]],unlist(obj[[split.by]]))
          
        }
      }
     else{
       if(!all(unique(unlist(obj[[split.by]]))%in%str_remove(Layers(obj[[allcount.assay]],search='counts'),'counts\\.'))){
         if(length(Layers(obj[[allcount.assay]],search='counts'))>1)
           obj[[allcount.assay]]<-JoinLayers(obj[[allcount.assay]])
         obj[[allcount.assay]]<-split(obj[[allcount.assay]],unlist(obj[[split.by]]))
         
     }
      
     }
      
      #calculate the splicing def matrix for each layer, saved in BPCells dir
      
      splicdef_list<-sapply(unique(unlist(obj[[split.by]])), function(ind){
        mat_dir=file.path(bpcells_dir,ind,assay.name)
        
        if(!file.exists(mat_dir)&!force_bpcells){
          message('Extracting unspliced and spliced matrices for ', ind)
          
          #get dgCmatrices
          count_layer=paste0('counts.',ind)
          
          
          us= LayerData(obj[[unspliced.assay]],count_layer)|>as(Class = 'dgCMatrix')
          
          if(is.null(spliced.assay)){
            s= LayerData(obj[[allcount.assay]],count_layer)|>as(Class = 'dgCMatrix')
            
          }else{
            s= LayerData(obj[[spliced.assay]],count_layer)|>as(Class = 'dgCMatrix')
            
          }
          
          #order,names and harmonize the matrices
          comm_genes=intersect(rownames(us),rownames(s))
          us=us[comm_genes,]
          s=s[comm_genes,]
          message('dimension of the unspliced matrix: ',paste(dim(us),collapse = 'x'))
          message('dimension of the spliced matrix: ',paste(dim(s),collapse = 'x'))
          if(any(dim(s)!=dim(us)))stop('different dimensions of the 2 matrices')
          #calc splic def
          message('Calculating splicing deficit..')
          
          if(is.null(spliced.assay)){
            splicingdef=CalcSplicDef(obj=us,all=s,...)
            
          }else{
            splicingdef=CalcSplicDef(obj=us,spliced=s,...)
            
          }
          
          message('saving Splicing deficit matrix in BPCells dir format..')
          
          
          if(dir.exists(mat_dir)){
            system(paste('rm -rf',mat_dir))
          }
          
          write_matrix_dir(
            mat =splicingdef,
            dir = mat_dir)
          message('done.')
          
          
        }
        
        splicingdef <- open_matrix_dir(dir =mat_dir)
        
        
        return(splicingdef)
        
      })
      
      #create the new assay
      message('Creating ',assay.name,' Assay')
      
      obj[[assay.name]] <- CreateAssay5Object(data = splicdef_list)
      
    }else{
      #Work with SeuratV3 splitting the object instead of creating layer
      #1) Split
      obj_list <- SplitObject(obj, split.by = split.by)
      inds=names(obj_list)
      obj_list=mclapply(obj_list,function(obj){
        ind=unique(as.vector(unlist(obj[[split.by]])))
        i=which(inds==ind)
        
        message('calculating splcing def for ',ind,' (',i,'/',length(inds),')')
        
        obj<-CalcSplicDef(obj,
                          unspliced.assay=unspliced.assay,
                          spliced.assay=spliced.assay,
                          allcount.assay=allcount.assay,
                          group.by=group.by,
                          split.by=NULL,on.disk=NULL,
                          ...)
        if(on.disk){
          mat_dir=file.path(bpcells_dir,ind,assay.name)
          if(dir.exists(mat_dir)){
            system(paste('rm -rf',mat_dir))
          }
          write_matrix_dir(
            mat = obj[[assay.name]]@data,
            dir = mat_dir)
          bpmat <- open_matrix_dir(dir =mat_dir)
          obj[[assay.name]] <- CreateAssay5Object(data = bpmat)
        }
        return(obj)
        
      },mc.cores=nThreads)
      
      #2) merge 
      message('merging the splitted objects')
      
      obj=merge(obj_list[[1]],obj_list[2:length(obj_list)])
      
      return(obj)
    }
    
    
  }else{
    #transform to Assay V3 object if not
    if('Assay5'%in%class(obj[[unspliced.assay]])){
      message('Assay 5 detected for unspliced assay, convert to Assay V3')
      
      # convert a v5 assay to a v3 assay
      obj[[unspliced.assay]] <- as(object = obj[[unspliced.assay]], Class = "Assay")
    }
    if(!is.null(spliced.assay)){
      
      if('Assay5'%in%class(obj[[spliced.assay]])){
        message('Assay 5 detected for spliced assay, convert to Assay V3')
        
        # convert a v5 assay to a v3 assay
        obj[[spliced.assay]] <- as(object = obj[[spliced.assay]], Class = "Assay")
      }
    }
    if(!is.null(allcount.assay)){
      
      if('Assay5'%in%class(obj[[allcount.assay]])){
        message('Assay 5 detected for ',allcount.assay,' assay, convert to Assay V3')
        
        # convert a v5 assay to a v3 assay
        obj[[allcount.assay]] <- as(object = obj[[allcount.assay]], Class = "Assay")
      }
    }
    
    #get the unspliced and count matrices
    message('Extracting unspliced and spliced matrices')
    
    us=obj[[unspliced.assay]]@counts
    if(is.null(spliced.assay)){
      
      s=obj[[allcount.assay]]@counts
      
    }else{
      s=obj[[spliced.assay]]@counts
      
    }
    
    #order,names and harmonize the matrices
    comm_genes=intersect(rownames(us),rownames(s))
    us=us[comm_genes,]
    s=s[comm_genes,]
    message('dimension of the unspliced matrix: ',paste(dim(us),collapse = 'x'))
    message('dimension of the spliced matrix: ',paste(dim(s),collapse = 'x'))
    if(any(dim(s)!=dim(us)))stop('different dimensions of the 2 matrices')
    #calc splic def
    if(is.null(spliced.assay)){
      splicingdef=CalcSplicDef(obj=us,all=s,...)
      
    }else{
      splicingdef=CalcSplicDef(obj=us,spliced=s,...)
      
    }
    #add as assay
    message('Creating ',assay.name,' Assay')
    obj[[assay.name]]<-CreateAssayObject(data=splicingdef)
    
  }
  
  #Finally, calc the splicing defect per cell group by cluster
  #- only cluster expressed-genes in >50% cells
  
  message('Calculating splicing deficit score')
  
  for(l in Layers(obj[[assay.name]],search ='data')){
    ind=str_remove(l,'data\\.')
    message('for ',ind)
    splicingdef=LayerData(obj[[assay.name]],l)|>as(Class = 'dgCMatrix')
    genes<-rownames(splicingdef)
    
    spliced=LayerData(obj[[ifelse(is.null(spliced.assay),allcount.assay,spliced.assay)]],str_replace(l,'data','counts'))|>as(Class = 'dgCMatrix')
    
    for(ct in unlist(unique(obj[[group.by]]))){
      if(!silent)message(ct,'...\n')
      cells<-intersect(colnames(splicingdef),colnames(obj)[unlist(obj[[group.by]]==ct)])
      message(length(cells),' cells')
      if(length(cells)>=min.cells){
        
        is.expr<-rowSums(spliced[genes,cells]>0)>0.5*length(cells)
        
        if(!silent)message(sum(is.expr),' genes are expressed (>50%)')
        
        genes_expr<-genes[is.expr]
        if(length(genes_expr)>=min.genes){
          obj@meta.data[cells,'avg.splicdef_group']<-apply(splicingdef[genes_expr,cells],2,function(x)mean(x[x!=Inf],na.rm=T))
          obj@meta.data[cells,'med.splicdef_group']<-apply(splicingdef[genes_expr,cells],2,function(x)median(x[x!=Inf],na.rm=T))
          
        }else{
          message('not passing min.genes threshold, returning NA for splicdef score')
          
        }
        
      }else{
        message('not passing min.cells threshold, returning NA for splicdef score')
        
      }
    }
    
  }
 
  return(obj)
}

#CalcSplicDef.matrix
#calculate the splicing ratio as a new matrix 
#obj: unspliced matrix
#spliced: spliced matrix
#all :  allcount (unspliced+spliced) matrix if spliced matrix not present
#exons_annot, introns_annot: path to ucsc extracted exons / introns regions (in bed format)
#transcripts_annot, gene_name_annot: path to conversion table between transcript_id and gene_id
# gene_name_annot: path to conversion table between gene_id and gene_name, needed if matrix rownames in gene_name format 


CalcSplicDef.Matrix<-function(obj,spliced=NULL,all=NULL,
                              exons_annot='/projectnb/tcwlab/RefData/ucsc/gencode.v44.knownGene.exons.clean.bed.gz',
                              introns_annot='/projectnb/tcwlab/RefData/ucsc/gencode.v44.knownGene.introns.clean.bed.gz',
                              transcripts_annot='/projectnb/tcwlab/RefData/ensembl/ensembl_transcript_gene_id_trans.csv.gz',
                              gene_name_annot='/projectnb/tcwlab/RefData/10XGenomics/features_gene_id_gene_name.tsv'){
  
  require(stringr)
  require(data.table)
  
  
  if(is.null(spliced)){
    message('calculating the spliced matrix..')
    spliced=all-obj
    message('done.')
  }
  #check if the matrices have the same size and order
  same.rn=all(rownames(spliced)==rownames(obj))
  same.cn=all(colnames(spliced)==colnames(obj))
  
  if(!(same.rn&same.cn))stop('the spliced and unspliced matrices do not have same row or columns order')
  
  have.gn=!all(str_detect(rownames(spliced),'^ENSG'))
  
  if(have.gn){
    message('detected gene names, converting as gene_id')
    transg<-fread(gene_name_annot)
    transg<-unique(transg,by='gene_name')
    #keep only genes present in the reference
    in.ref=rownames(obj)%in%transg$gene_name
    message(sum(in.ref),' Genes found matching reference, convert them')
    obj<-obj[in.ref,]
    spliced<-spliced[in.ref,]
    
    #then translate to gene_id
    rownames(obj)<-transg[rownames(obj),on='gene_name']$gene_id
    rownames(spliced)<-transg[rownames(spliced),on='gene_name']$gene_id
    
    message('translated rownames: ',paste(head(rownames(obj)),collapse = ', '))
  }
  
  
  #keep genes with exons and introns detected in >5% of the cells 
  exons_expressed<-rowSums(spliced>0,)>0.05*ncol(spliced)
  introns_expressed<-rowSums(obj>0)>0.05*ncol(obj)
  
  genes_id_to_keep<-intersect(rownames(spliced)[exons_expressed],
                              rownames(obj)[introns_expressed])
  
  message(length(genes_id_to_keep),' genes pass 5% expression QC threshold')
  
  #normalize intron (unspl) and exons (splice) count by their size
  trans<-fread(transcripts_annot)
  setnames(trans,colnames(trans),str_remove(colnames(trans),'^ensembl_'))
  
  #get introns/exons size
  anno_exons<-fread(exons_annot,col.names = c('chr','start','end','ID'))
  anno_exons[,size:=end-start]
  anno_exons[,transcript_id:=str_extract(ID,'ENST[0-9]+')]
  anno_exons[,size_by_transcript:=sum(size),by='transcript_id']
  
  
  anno_exons<-merge(anno_exons,trans)
  
  
  #int size
  anno_introns<-fread('~/tcwlab/RefData/ucsc/gencode.v44.knownGene.introns.clean.bed.gz',col.names = c('chr','start','end','ID'))
  anno_introns[,size:=end-start]
  anno_introns[,transcript_id:=str_extract(ID,'ENST[0-9]+')]
  anno_introns[,size_by_transcript:=sum(size),by='transcript_id']
  
  anno_introns<-merge(anno_introns,trans)
  
  #keep one transcript per gene 
  anno_exonsf<-unique(anno_exons,by='gene_id')
  
  anno_intronsf<-unique(anno_introns,by='gene_id')
  

  genes_id_to_keep<-intersect(genes_id_to_keep,anno_intronsf$gene_id)
  genes_id_to_keep<-intersect(genes_id_to_keep,anno_exonsf$gene_id)
  
  
  message('Calculating the unspliced ratio matrix for ',length(genes_id_to_keep),' genes')
  
  #calc splicdef =(count introns/size introns)/(count exons/size exons)
  
  splicingdef<-(obj[genes_id_to_keep,]/anno_intronsf[genes_id_to_keep,on='gene_id']$size_by_transcript)/
    (spliced[genes_id_to_keep,]/anno_exonsf[genes_id_to_keep,on='gene_id']$size_by_transcript)
  rownames(splicingdef)<-genes_id_to_keep
  colnames(splicingdef)<-colnames(obj)
  
  #transform as dgCMatrix
  splicingdef<-as(splicingdef,Class = 'dgCMatrix')
  
  if(have.gn){
    message('trasnform back to gene names')
    rownames(splicingdef)<-transg[rownames(splicingdef),on='gene_id']$gene_name
  }
  
  return(splicingdef)
  
  
}
