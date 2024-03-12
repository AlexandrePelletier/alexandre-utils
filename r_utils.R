#R utils packages and functions for every project

#constantly used packages####
library("data.table")
setDTthreads(threads = 0)
message(getDTthreads(),' threads available for data.table')

library("stringr")
library("ggplot2")
library("ggrepel")
library("patchwork")
library("here")


#basic utils functions####

fp<-function(...)file.path(...)

ps<-function(...,sep="",collapse = NULL)paste(...,sep=sep,collapse = collapse)

source<-function(file,chdir=TRUE)base::source(file,chdir = chdir)

# Stats Related Functions ####
CategoricalsToDummy<-function(covs){
  covs2<-copy(covs)
  id_col<-colnames(covs2)[as.vector(unlist(covs2[,lapply(.SD, function(x)length(unique(x))==.N)]))]
  num_vars<-colnames(covs2[,.SD,.SDcols=is.numeric])
  cat_vars<-setdiff(colnames(covs2)[unlist(covs2[,lapply(.SD, function(x)all(table(x)>1))])],num_vars)
  others_vars<-setdiff(colnames(covs2),c(num_vars,cat_vars,id_col))
  message(paste(others_vars,collapse = ', '),' are categorical variables with only 1 sample falling in one category, so wont be Dummyfied')
  dum_vars<-c()
  
  for(cat in cat_vars){
    t<-table(as.factor(unlist(covs2[,..cat])))
    lvls<-names(t)[-1]
    cols<-paste0(cat,lvls)
    covs2[,(cols):=lapply(lvls, function(l)as.numeric(.SD==l)),.SDcols=cat]
    
    dum_vars<-c(dum_vars,cols) 
  }
  
  return(covs2[,.SD,.SDcols=c(id_col,num_vars,others_vars,dum_vars)])
}



GetVarPCs<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}


#Reformatting classical R data/results ####
UniqueClean<-function(x,key_cols='sample_id',pattern_to_exclude=NULL){
  #keep only unconstant variables
  nums_to_keep<-names(which(x[,sapply(.SD,function(y)(!all(is.na(y)))&var(y,na.rm = T)!=0)&length(unique(y))!=.N,.SDcols=is.numeric]))
  
  cats_to_keep<-names(which(x[,sapply(.SD,function(y)(!all(is.na(y)))&length(unique(y))!=.N&length(unique(y))!=1),.SDcols=!is.numeric]))
  
  if(!is.null(pattern_to_exclude)){
    nums_to_keep<-nums_to_keep[!str_detect(nums_to_keep,pattern_to_exclude)]
    cats_to_keep<-cats_to_keep[!str_detect(cats_to_keep,pattern_to_exclude)]
    
  }
  cols_to_keep<-c(key_cols,nums_to_keep,cats_to_keep)
  
  return(unique(x,by=c(key_cols))[,.SD,.SDcols=cols_to_keep])
}

RemoveUselessColumns<-function(x,key_cols='sample_id',pattern_to_exclude=NULL){
  #keep only unconstant variables
  nums_to_keep<-names(which(x[,sapply(.SD,function(x)(!all(is.na(x)))&var(x,na.rm = T)!=0),.SDcols=is.numeric]))
  
  cats_to_keep<-names(which(x[,sapply(.SD,function(x)(!all(is.na(x)))&length(unique(x))!=.N&length(unique(x))!=1),.SDcols=!is.numeric]))
  
  if(!is.null(pattern_to_exclude)){
    nums_to_keep<-nums_to_keep[!str_detect(nums_to_keep,pattern_to_exclude)]
    cats_to_keep<-cats_to_keep[!str_detect(cats_to_keep,pattern_to_exclude)]
    
  }
  cols_to_keep<-c(key_cols,nums_to_keep,cats_to_keep)
  
  return(x[,.SD,.SDcols=cols_to_keep])
}

DetectDEResFormat<-function(res_de){
  
  if(all(c('log2FoldChange','padj')%in%colnames(res_de)))return('DESEQ2')
  else if(any(c('avg_log2FC','avg_logFC')%in%colnames(res_de))&'p_val_adj'%in%colnames(res_de))return('SEURAT')
  else return('unknown')
  
}

FormatDEResToSeurat<-function(res_de){
  res_forms<-list(DESEQ2=c(padj='padj',
                           pval='pvalue',
                           stat='stat',
                           FC='log2FoldChange'),
                  SEURAT=c(padj='p_val_adj',
                           pval='p_val',
                           stat='stat',
                           FC='avg_log2FC'))
  
  if(DetectDEResFormat(res_de)=='DESEQ2'){
    setnames(res_de,old = res_forms[['DESEQ2']],
             new = res_forms[['SEURAT']])
    return(res_de)
  }else{
    stop('unknown format')
  }
  return(res_de)
}

ReFormatDERes<-function(res_de,to='SEURAT'){
  if(to=='SEURAT') ifelse(DetectDEResFormat(res_de)==to,return(res_de),return(FormatDEResToSeurat(res_de)))
  else{
    stop('unsupported formatting')
  }
}


#GGPLOT####

bar_bw<-function()scale_fill_manual(values=c('black','grey'))
bar_rb<-function(invert=FALSE)ifelse(invert,return(scale_fill_manual(values=c('royalblue3','orangered3'))),return(scale_fill_manual(values=c('orangered3','royalblue3'))))

#BIOMART####


GetMartGenes<-function()biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
GetMartMouseGenes<-function()biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

GetMartReg<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")
GetMartMotif<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_motif_feature")



GetBiomartAttrs<-function(mart)data.table::data.table(biomaRt::listAttributes(mart))
GetBiomartFilter<-function()data.table::data.table(biomaRt::listFilter(mart))

GetBiomartAttrs_reg<-function()data.table::data.table(listAttributes(biomaRt::useEnsembl(biomart = "regulation",dataset = "hsapiens_regulatory_feature")))
GetBiomartFilter_reg<-function()data.table::data.table(biomaRt::listFilter(biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")))

TransNMtoSymbol<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
                                      filters = 'refseq_mrna', 
                                      values = refseq_ids, 
                                      mart = GetMartGenes())))
}

TransTranscriptToGene<-function(transcript_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'),
                                               filters = 'ensembl_transcript_id', 
                                               values = transcript_ids, 
                                               mart = GetMartGenes())))
}

TransSymboltoNM<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}

TransEnsembltoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                                               filters = 'ensembl_gene_id', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransSymboltoEnsembl<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}

TransEnsemblVerstoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
                                               filters = 'ensembl_gene_id_version', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransSymboltoEnsemblVers<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}



TransEnsembltoNM<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna'),
                                               filters = 'ensembl_gene_id', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransNMtoEnsembl<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna'),
                                               filters = 'refseq_mrna', 
                                               values = refseq_ids, 
                                               mart = GetMartGenes())))
}

TransEnsemblVerstoNM<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna'),
                                               filters = 'ensembl_gene_id_version', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransNMtoEnsemblVers<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna'),
                                               filters = 'refseq_mrna', 
                                               values = refseq_ids, 
                                               mart = GetMartGenes())))
}



tr<-function(ids_sepBySlash,retourne="all",sep="/",tradEntrezInSymbol=FALSE,uniqu=TRUE){
  IDs<-as.vector(strsplit(ids_sepBySlash,sep)[[1]])
  if(retourne=="all"){
    ret<-1:length(IDs)
  }else{
    ret<-retourne
  }
  if(tradEntrezInSymbol){
    require(clusterProfiler)
    library(org.Hs.eg.db)
    if(retourne=="all"){
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL)
    }
    else{
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL[ret])
    }
    
  }else{
    if(uniqu){
      return(unique(IDs[ret]))
    }else{
      return(IDs[ret])
    }
    
    
  }
  
}







# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x,return_dt=T){
  
  require("biomaRt")
  human =GetMartGenes()
  mouse = GetMartMouseGenes()
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  if(return_dt)return(data.table(genesV2))
  
  mousex <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

convertMouseGeneList <- function(x,return_dt=T){
  
  require("biomaRt")
  human =GetMartGenes()
  mouse = GetMartMouseGenes()
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  if(return_dt)return(data.table(genesV2))
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


#Over repre / GSEA####


OR<-function(set1,set2,size_universe){
  if(any(duplicated(set1))){
    set1<-unique(set1)
  }
  if(any(duplicated(set2))){
    set2<-unique(set2)
  }
  phyper(q=sum(set1%in%set2)-1, 
         #number of white balls drawn without replacement from an urn which contains both black and white balls. 
         #<=> number of tcf4 target in DEGs. "-1" because normally give P[X > x] but we want P[X >= x])
         m=length(set2), #the number of white balls in the urn. <=> number of DEGs 
         n=size_universe-length(set2), #the number of black balls in the urn. <=> number of genes tested - number of DEGs
         k=length(set1), #the number of balls drawn from the urn <=> numbers of tcf4 targets 
         lower.tail=FALSE)  # if TRUE (default), probabilities are P[X ≤ x] (under-representation), otherwise, P[X > x] (over-representation).
  
}



OR2<-function(querys,terms_list,size_universe,min.term.size=0,max.term.size=Inf,overlap_column=TRUE,verbose=FALSE){
  if(is.list(querys)){
    return(Reduce(rbind,lapply(names(querys),
                               function(q)OR2(querys = querys[[q]],
                                              terms_list = terms_list,
                                              size_universe = size_universe,
                                              min.term.size = min.term.size,
                                              max.term.size = max.term.size,
                                              overlap_column = overlap_column,
                                              verbose = verbose )[,query:=q])))
  }else{
    res_or<-data.table(term=names(terms_list),term.size=sapply(terms_list,length))
    res_or<-res_or[term.size<=max.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(length(terms_list)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
    res_or<-res_or[term.size>=min.term.size]
    if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
    
    res_or[,n.query:=length(querys)]
    res_or[,n.overlap:=sum(querys%in%terms_list[[term]]),by="term"]
    if(overlap_column==TRUE){
      res_or[,genes.overlap:=paste(querys[querys%in%terms_list[[term]]],collapse="|"),by="term"]
    }
    res_or[,pct.query.overlap:=n.overlap/n.query]
    res_or[,precision:=pct.query.overlap]
    
    res_or[,pct.term.overlap:=n.overlap/term.size]
    
    res_or[,background_size:=size_universe]
    
    res_or[,pct.background:=term.size/size_universe] #TO IMPROVE (Here size universe can be the intersection between 2 "univers" while n.query is the n.query is the query universe. 
    
    
    res_or[,pval:=phyper(q=n.overlap-1, 
                         m=n.query, 
                         n=size_universe-n.query, 
                         k=term.size, 
                         lower.tail=FALSE),
           by="term"]
    res_or[,padj:=p.adjust(pval,method = 'BH')]
    if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
    return(res_or)
  }
  
}

OR3<-function(querys,terms_list,background,min.term.size=0,max.term.size=Inf,overlap_column=TRUE,verbose=FALSE){
  if(is.list(querys)){
    dt<-Reduce(rbind,lapply(names(querys),
                            function(q)OR3(querys = querys[[q]],
                                           terms_list = terms_list,
                                           background = background,
                                           min.term.size = min.term.size,
                                           max.term.size = max.term.size,
                                           overlap_column = overlap_column,
                                           verbose = verbose )[,query:=q]))
    
    return(dt[,query.:=query][,.SD,.SDcols=c(ncol(dt),1:(ncol(dt)-1))])
  }else{
    queryf<-intersect(querys,background)
    terms_listf<-lapply(terms_list,function(x)intersect(x,background))
    res_or<-data.table(term=names(terms_listf),term.size=sapply(terms_listf,length))
    res_or<-res_or[term.size<=max.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(length(terms_listf)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
    res_or<-res_or[term.size>=min.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
    
    res_or[,n.query:=length(queryf)]
    res_or[,n.overlap:=sum(queryf%in%terms_listf[[term]]),by="term"]
    
    res_or[,pct.query.overlap:=n.overlap/n.query]
    res_or[,precision:=pct.query.overlap]
    
    res_or[,pct.term.overlap:=n.overlap/term.size]
    
    res_or[,background_size:=length(background)]
    
    res_or[,pct.term.background:=term.size/background_size] 
    
    
    res_or[,pval:=phyper(q=n.overlap-1, 
                         m=term.size, 
                         n=background_size-term.size, 
                         k=n.query, 
                         lower.tail=FALSE),
           by="term"]
    res_or[,padj:=p.adjust(pval,method = 'BH')]
    res_or[,fold.enrichment:=pct.query.overlap/pct.term.background]
    if(overlap_column==TRUE){
      res_or[,genes.overlap:=paste(queryf[queryf%in%terms_listf[[term]]],collapse="|"),by="term"]
    }
    if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
    return(res_or)
  }
  
}

#trans in hg38
hg19to38<-function(x){
  in_file<-"outputs/temp_hg19.bed"
  out_file<-"outputs/temp_hg38.bed"
  
  fwrite(x,in_file,col.names = F,sep="\t")
  
  system(paste("CrossMap.py bed ref/hg19ToHg38.over.chain.gz",in_file,out_file))
  trans<-fread(out_file,select=c(1,2,3,4),col.names = c("chr","start","end","id"))
  file.remove(c(in_file,out_file))
  return(trans)
}


FindGOGenes<-function(terms_or_ids){
  require("biomaRt")
  require("stringr")
  
  if(!all(str_detect(terms_or_ids,"^GO:")))terms_or_ids=FindGO_ID(term_description=terms_or_ids)
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
  gene.data <- getBM(attributes=c('hgnc_symbol','go_id'),
                     filters = 'go', values = terms_or_ids, mart = ensembl,uniqueRows = T)
  return(data.table(gene.data)[go_id%in%terms_or_ids])
}



FindGO_ID<-function(term_descriptions){
  require("GO.db")
  terms<-Term(GOTERM)
  ids<-names(terms)[match(term_descriptions,unlist(terms))]
  return(ids)
}




####Signac####
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

CheckMotif<-function(object,peaks,motif.name,assay = NULL,return.peaks=FALSE){
  require("Signac")
  if(is.null(assay))assay<-DefaultAssay(object)
  motif<-GetMotifIDs(object,motif.name,assay=assay)
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  
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
GetMotif<-function(object,peaks,motifs=NULL,assay = "peaks"){
  #return tfmotif peak data.table
  require("Signac")
  require("data.table")
  
  if(is.null(assay))assay<-DefaultAssay(object)
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
  
  motifsnames<-data.table(motif.name=object@assays$peaks@motifs@motif.names,motif=names(object@assays$peaks@motifs@motif.names))
  motif_dt<-merge(motif_dt,motifsnames)
  return(motif_dt)
}




start<-function(x)sapply(x,function(x)as.numeric(strsplit(x,"-")[[1]][2]))
end<-function(x)sapply(x,function(x)as.numeric(strsplit(x,"-")[[1]][3]))
seqid<-function(x)sapply(x,function(x)strsplit(x,"-")[[1]][1])


MethChangeReg<-function(res_meth,region){
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
  res_meth_reg[,start:=pos][,end:=pos+1]
  return(res_meth_reg)
}
MethChangePlot<-function(res_meth,region,limits=NULL,breaks=waiver()){
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
  res_meth_reg[,start:=pos][,end:=pos+1]
  p<-ggplot(data = res_meth_reg) + geom_segment(aes(x = start, y = 0, 
                                                    xend = end, yend = logFC,col=-log10(P.Value)), size = 2, data = res_meth_reg)+
    scale_color_gradient(low = "white",high = "black",limits=limits,breaks=breaks)
  
  p<-p+ theme_classic() + ylab(label = "Methylation change") + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  
  return(p)
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

TFsMotifPlot<-function(object,region,motif.names,assay=NULL,size=2,alpha=1,pad=0){
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
  
  
  p<-ggplot(data = dt_region) + geom_segment(aes(x = start-pad, y = 0, 
                                                 xend = end+pad, yend = 0,col=motif.name), size = size,alpha=alpha, data = dt_region)
  
  p<-p+ theme_classic() + ylab(label = "TF motif") + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  return(p)
}



str_to_vec<-function(ids_sepBySlash,sep="/")as.vector(strsplit(ids_sepBySlash,sep)[[1]])

#DATA.TABLE ADD IN
freadvcf<-function(file){
  cmd<-paste('zcat',file,'| grep -v "^##"')
  fread(cmd)
}


#GSUB FILES CREATION####
CreateJobFile<-function(cmd_list,file,proj_name='tcwlab',modules=NULL,
                        loadBashrc=FALSE,conda_env=NULL,
                        micromamba_env=NULL,
                        nThreads=NULL,memPerCore=NULL,maxHours=24,
                        parallelize=FALSE){
  
  template_header='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_header.txt'
  template_tail='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_tail.txt'
  filename<-basename(file)
  projdir<-dirname(file)
  
  while(str_detect(projdir,'scripts')){
    projdir<-dirname(projdir)
  }
  
  dir.create(file.path(projdir,'logs'),showWarnings = F)
  file_path<-file
  log_file=file.path(projdir,'logs',paste0(str_remove(filename,'.qsub$'),'.log'))
  
  #create the qsub file
  cat('#!/bin/bash -l\n',file = file_path)
  message('qsub file created at ',file_path)
  
  #set the general job parameters
  proj_name_opt<-paste('-P ',proj_name)
  CombStdOutErr_opt<-'-j y'
  maxHours_opt<-paste0('-l h_rt=',maxHours,':00:00')
  qlog<-paste('-o ',log_file)
  if(!is.null(nThreads)){
    nThreads_opt=paste('-pe omp',nThreads)
  }else{
    nThreads_opt=NULL
    
  }
  
  if(!is.null(memPerCore)){
    memPerCore_opt=paste0('-l mem_per_core=',memPerCore)
  }else{
    memPerCore_opt=NULL
    
    }
  #create child qsub files if parallelize mode
  
  if(parallelize&length(cmd_list)>1){
    
    script_id=str_extract(filename,'^[0-9A-Za-z]+')
    scripts_dir<-file.path(projdir,'scripts')
   if(!dir.exists(scripts_dir))dir.create(scripts_dir)
    
    #add the main job parameters
    cat( '#Parameters of the Jobs :',file = file_path,append = T)
    cat( c('\n',proj_name_opt,CombStdOutErr_opt,maxHours_opt,qlog),
         file = file_path,append = T,sep = '\n#$')
    cat( '\n',file = file_path,append = T,sep = '\n')
    
    
      if(!is.null(names(cmd_list))){
        child_jobnames<-make.names(names(cmd_list))
      }else{
        child_jobnames<-paste0('child',1:length(cmd_list))
      }
    
    child_jobfiles<-file.path(scripts_dir,paste0(script_id,'child-',child_jobnames,'.qsub'))
   
     #add command to run the childs qsub
    cmds<-sapply(1:length(cmd_list),function(i){
      
      cmd<-cmd_list[[i]]
      CreateJobFile(cmd_list =cmd,file = child_jobfiles[[i]],nThreads = nThreads ,proj_name = proj_name ,
                    loadBashrc = loadBashrc,modules = modules,conda_env = conda_env,micromamba_env = micromamba_env,
                    maxHours =  maxHours,memPerCore = memPerCore,parallelize = FALSE)
      
      # Submit the job using qsub and capture the job ID
      cmd<-paste(paste0('job_id',i),'=$(qsub', '-N',paste0('j',script_id,child_jobnames[[i]]),child_jobfiles[[i]],proj_name)
      # cmd<-RunQsub(child_jobfiles[[i]],job_name = paste0('j',script_id,child_jobnames[[i]]),dryrun = T)
      return(cmd)
      
    })
    
    cat( cmds,file = file_path,append = T,sep = '\n')
    
    # Wait until the jobs are completed
    for(i in 1:length(cmd_list)){
      cmd_wait<-c(paste('while qstat | grep -w',paste0("$job_id",i),  '> /dev/null; do'),
                  paste('echo', '"waiting jobs to complete.."'),
                  paste('sleep', '30'), # Adjust sleep time as needed
                  'done',
                  paste('echo',paste0('"Job',i),'have completed."'))
      cat( c('\n',cmd_wait,file = file_path,append = T,sep = '\n'))
      
    }
   
    
  }else{
    #add the core job parameters
  cat( '#Parameters of the Jobs :',file = file_path,append = T)
  cat( c('\n',proj_name_opt,CombStdOutErr_opt,maxHours_opt,qlog,nThreads_opt,memPerCore_opt),file = file_path,append = T,sep = '\n#$')
  cat( '\n',file = file_path,append = T,sep = '\n')
  
  #add the module to loads
  if(!is.null(modules)){
    modules<-ifelse(modules=='R','R/4.2.1',modules)
    
    if('gatk'%in%modules){
      
      modules<-c(modules[modules!='gatk'],
                 'miniconda/23.1.0',
            'java/17.0.8',
            'gatk/4.4.0.0')
      
      conda_env<-union(conda_env,'/share/pkg.8/gatk/4.4.0.0/install/gatk-4.4.0.0')
      
    }
    
    cat( '#Modules to load:',file = file_path,append = T)
    cat( c('\n',modules),file = file_path,append = T,sep = '\nmodule load ')
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  #load .bashrc
  if(loadBashrc){
    cat( '# loading of bashrc profile:',file = file_path,append = T)
    cat('\n source $HOME/.bashrc',file = file_path,append = T)
    cat( '\n',file = file_path,append = T,sep = '\n')
  }
  
  #activate conda environment 
  if(!is.null(conda_env)){
    cat( '#Conda environment activation:',file = file_path,append = T)
    cat( c('\n',conda_env),file = file_path,append = T,sep = '\nconda activate ')
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  
  #activate micromamba environment 
  if(!is.null(micromamba_env)){
    cat( '#Micromamba environment activation:',file = file_path,append = T)
    cat( c('\n',micromamba_env),file = file_path,append = T,sep = '\nmicromamba activate ')
    if('pisces-rabbit'%in%micromamba_env){
      cat( c('\n# singularity specifics configuration: ',
             'export SINGULARITY_BIND="$TMP,/restricted/projectnb/tcwlab-adsp/,/projectnb/tcwlab-adsp/,/projectnb/tcwlab/"',
             'touch .Rprofile'),
           file = file_path,append = T,sep = '\n')

    }
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  
  
  #add the header
  system(paste('cat',template_header,'>>',file_path))
  
  #add the commandes to exec

  if(!is.null(names(cmd_list))){
    cmds<-unlist(lapply(names(cmd_list),
                        function(s)return(c(paste('echo',paste0('"----- Processing of ',s,' -----"')),cmd_list[[s]]))))
    
  }else{
    cmds<-unlist(cmd_list)
  }
  
  cat( cmds,file = file_path,append = T,sep = '\n')
  
  if('pisces-rabbit'%in%micromamba_env){
    cat( c('\n# Singularity specifics clean up: ',
           'rm .Rprofile'),file = file_path,append = T,sep = '\n')
  }
  cat( '\n',file = file_path,append = T,sep = '\n')
  
  #add the tail
  system(paste('cat',template_tail,'>>',file_path))
  
  #show the file header
  message('header:')
  system(paste('head -n 15',file_path))
  }
  
  #show the first commands
  message('5 first commands:')
  cat( head(cmds,5),sep = '\n')
  
}

CreateJobForPyFile<-function(python_file,proj_name='tcwlab',modules=NULL,
                        loadBashrc=FALSE,conda_env=NULL,
                        micromamba_env=NULL,
                        nThreads=4,memPerCore=NULL,maxHours=24){
  template_header='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_header.txt'
  template_tail='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_tail.txt'
  
  filename<-basename(python_file)
  projdir<-ifelse(str_detect(python_file,'scripts/'),dirname(dirname(python_file)),dirname(python_file))
  
  
  #create qsub file
  qsub_file<-str_replace(python_file,'\\.py$','.qsub')
  cat('#!/bin/bash -l\n',file = qsub_file)
  
  #create log file
  dir.create(file.path(projdir,'logs'),showWarnings = F)
  log_file=file.path(projdir,'logs',str_replace(filename,'\\.py$','.log'))
  
  
  #add the job parameters
  proj_name<-paste('-P ',proj_name)
  CombStdOutErr<-'-j y'
  maxHours<-paste0('-l h_rt=',maxHours,':00:00')
  qlog<-paste('-o ',log_file)
  if(!is.null(nThreads))nThreads=paste('-pe omp',nThreads)
  if(!is.null(memPerCore))memPerCore=paste0('-l mem_per_core=',memPerCore)
  
  cat( '#Parameters of the Jobs :',file = qsub_file,append = T)
  cat( c('\n',proj_name,CombStdOutErr,maxHours,qlog,nThreads,memPerCore),file = qsub_file,append = T,sep = '\n#$')
  cat( '\n',file = qsub_file,append = T,sep = '\n')
  
  #add the module to loads
  if(!is.null(modules)){
    modules<-ifelse(modules=='R','R/4.2.1',modules)
    
    if('gatk'%in%modules){
      
      modules<-c(modules[modules!='gatk'],
                 'miniconda/23.1.0',
                 'java/17.0.8',
                 'gatk/4.4.0.0')
      
      conda_env<-union(conda_env,'/share/pkg.8/gatk/4.4.0.0/install/gatk-4.4.0.0')
      
    }
    
    cat( '#Modules to load:',file = qsub_file,append = T)
    cat( c('\n',modules),file = qsub_file,append = T,sep = '\nmodule load ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  #load .bashrc
  if(loadBashrc){
    cat( '# loading of bashrc profile:',file = qsub_file,append = T)
    cat('\n source $HOME/.bashrc',file = qsub_file,append = T)
    cat( '\n',file = qsub_file,append = T,sep = '\n')
  }
  
  #activate conda environment 
  if(!is.null(conda_env)){
    cat( '#Conda environment activation:',file = qsub_file,append = T)
    cat( c('\n',conda_env),file = qsub_file,append = T,sep = '\nconda activate ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  
  #activate micromamba environment 
  if(!is.null(micromamba_env)){
    cat( '#Micromamba environment activation:',file = qsub_file,append = T)
    cat( c('\n',micromamba_env),file = qsub_file,append = T,sep = '\nmicromamba activate ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  

  #add the header
  system(paste('cat',template_header,'>>',qsub_file))
  
  
  #add the Python script to exec
  cmd<-paste('python',python_file,'>>',log_file)
  
  cat( cmd,file = qsub_file,append = T,sep = '\n')
  cat( '\n',file = qsub_file,append = T,sep = '\n')
  
  #add the tail
  system(paste('cat',template_tail,'>>',qsub_file))
  
  #show the file header
  message('qsub file created at ',qsub_file)
  message('header:')
  system(paste('head -n 15',qsub_file))
  
  #show the  bash command
  message('the bash command to execute:')
  cat( cmd,sep = '\n')
  
  #show the first R lines
  message('the 15 firsts lines of Python scripts to execute:')
  system(paste('head -n 15',python_file))
  
  
  
}


CreateJobForRfile<-function(r_file,proj_name='tcwlab',
                            modules='R',
                            loadBashrc=FALSE,conda_env=NULL,micromamba_env=NULL,
                            nThreads=4,memPerCore=NULL,maxHours=24){
  
  template_header='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_header.txt'
  template_tail='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_tail.txt'
  
  filename<-basename(r_file)
  projdir<-ifelse(str_detect(r_file,'scripts/'),dirname(dirname(r_file)),dirname(r_file))
  
  
  #create qsub file
  qsub_file<-str_replace(r_file,'\\.R$','.qsub')
  cat('#!/bin/bash -l\n',file = qsub_file)
  
  #create log file
  dir.create(file.path(projdir,'logs'),showWarnings = F)
  log_file=file.path(projdir,'logs',str_replace(filename,'\\.R$','.log'))

  
  #add the job parameters
  proj_name<-paste('-P ',proj_name)
  CombStdOutErr<-'-j y'
  maxHours<-paste0('-l h_rt=',maxHours,':00:00')
  qlog<-paste('-o ',log_file)
  if(!is.null(nThreads))nThreads=paste('-pe omp',nThreads)
  if(!is.null(memPerCore))memPerCore=paste0('-l mem_per_core=',memPerCore)
  
  cat( '#Parameters of the Jobs :',file = qsub_file,append = T)
  cat( c('\n',proj_name,CombStdOutErr,maxHours,qlog,nThreads),file = qsub_file,append = T,sep = '\n#$')
  cat( '\n',file = qsub_file,append = T,sep = '\n')
  
  #add the module to loads
  if(!is.null(modules)){
    modules<-ifelse(modules=='R','R/4.2.1',modules)
    cat( '#Modules to load:',file = qsub_file,append = T)
    cat( c('\n',modules),file = qsub_file,append = T,sep = '\nmodule load ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  
  #load .bashrc
  if(loadBashrc){
    cat('\n source $HOME/.bashrc',file = qsub_file,append = T)
    cat( '\n',file = qsub_file,append = T,sep = '\n')
  }
  
  #activate conda environment 
  if(!is.null(conda_env)){
    cat( '#Conda environment activation:',file = qsub_file,append = T)
    cat( c('\n',conda_env),file = qsub_file,append = T,sep = '\nconda activate ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  
  #activate micromamba environment 
  if(!is.null(micromamba_env)){
    cat( '#Micromamba environment activation:',file = qsub_file,append = T)
    cat( c('\n',micromamba_env),file = qsub_file,append = T,sep = '\nmicromamba activate ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  
  
  
  #add the header
  system(paste('cat',template_header,'>>',qsub_file))
  
  
  #add the Rscript to exec
  cmd<-paste('Rscript',r_file,'>>',log_file)
  
  cat( cmd,file = qsub_file,append = T,sep = '\n')
  cat( '\n',file = qsub_file,append = T,sep = '\n')
  
  #add the tail
  system(paste('cat',template_tail,'>>',qsub_file))
  
  #show the file header
  message('qsub file created at ',qsub_file)
  message('header:')
  system(paste('head -n 15',qsub_file))
  
  #show the  bash command
  message('the bash command to execute:')
  cat( cmd,sep = '\n')
  
  #show the first R lines
  message('the 15 firsts R lines to executes')
  system(paste('head -n 15',r_file))
  
  
}

RunQsub<-function(qsub_file,job_name,proj_name=NULL,wait_for=NULL,dryrun=FALSE){
  if(!str_detect(qsub_file,'.qsub$'))qsub_file=paste0(tools::file_path_sans_ext(qsub_file),'.qsub')
  
  if(!file.exists(qsub_file))stop(qsub_file,' do not exist')
  
  if(is.null(wait_for)){
    cmd<-paste('qsub','-N',job_name,qsub_file,proj_name)
 
  }else{
    cmd<-paste('qsub',
               '-N',job_name,
               '-hold_jid',wait_for,
               qsub_file,proj_name)
  }
  if(!dryrun){
    message<-system(cmd,intern = TRUE)
    message(paste(message,collapse = '\n'))
    
    jobid<-ifelse(str_detect(message,'has been submitted'),str_extract(message,'[0-9]+'),NA)
    return(jobid[!is.na(jobid)])
  }else{
    return(cmd)
  }

}



WaitQsub<-function(qsub_file,jobid,max_hours=24){
  if(!str_detect(qsub_file,'\\.qsub$'))qsub_file=paste0(tools::file_path_sans_ext(qsub_file),'.qsub')
  
  if(!file.exists(qsub_file))stop(qsub_file,' file does not exist')
  
  filename<-basename(qsub_file)
  projdir<-dirname(qsub_file)
  t0 <- Sys.time()
  
  while(str_detect(projdir,'scripts')){
    projdir<-dirname(projdir)
  }
  log_file=file.path(projdir,'logs',paste0(str_remove(filename,'.qsub$'),'.log'))
  
  max_sec=max_hours*60*60
  
  
  while(!file.exists(log_file)){
    Sys.sleep(60)
    t1 <- Sys.time()
    d<-t1-t0
    if(d>max_sec)stop('time reached max limits')
    message('waiting job ',jobid,' to be executed..')
    
    
  }
  
  jobid_inlog=str_extract(tail(grep('job ID :',
                                    readLines(con = log_file,skipNul = TRUE),
                                    value = T),n = 1),'[0-9]+$')

  while(jobid!=jobid_inlog){
    message('waiting job ',jobid,' to be executed..')
    Sys.sleep(60)
    t1 <- Sys.time()
    d<-t1-t0
    if(d>max_sec)stop('time reached max limits')
    jobid_inlog=str_extract(tail(grep('job ID :',
                                      readLines(con = log_file,skipNul = TRUE),
                                      value = T),n = 1),'[0-9]+$')
    
    
  }
  
  pattern=paste('Finished Analysis for job',jobid)
  lastlines<-tail(readLines(con = log_file,skipNul = TRUE))
  
  while(!any(str_detect(lastlines,pattern)) ){
    message('waiting job ',jobid,' to be completed..')
    Sys.sleep(120)
    lastlines<-tail(readLines(con = log_file,skipNul = TRUE))
    t1 <- Sys.time()
    d<-t1-t0
    if(d>max_sec)stop('time reached max limits')
  }
  t1 <- Sys.time()
  d<-t1-t0
  message(pattern,' after ', round(d/360,2),' hours')
  
  
}


#bash tools wrapper####
bed_inter<- function(a, b, opt1="-wa", opt2="-wb",out_dir=".", select=NULL, col.names=NULL){
  require(data.table)
  l<-list(a,b)
  files_to_rm<-c(FALSE,FALSE)
  file_paths<-sapply(1:2, function(i){
    x<-l[[i]]
    if(is.data.frame(x)){
      files_to_rm[i]<-TRUE
      file_path<-fp(out_dir,paste0("temp",i,".bed"))
      fwrite(x,file_path,sep="\t",col.names = FALSE)
      
    }else{
      file_path<-x
    }
    return(file_path)
  })
  
  out_file<-fp(out_dir,"temp_inter.bed")
  cmd<-paste("bedtools intersect -a",file_paths[1],"-b",file_paths[2], opt1, opt2,">",out_file)
  message("run in shell : ",cmd)
  system(cmd)
  message("done.")
  
  if(!is.null(col.names)){
    dt<-fread(out_file,select = select,col.names = col.names)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }else{
    dt<-fread(out_file,select = select)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }
}


RunGatk<-function(cmd){
  message('run in a terminal:')
  
  cat(paste('module load miniconda/23.1.0',
            'module load java/17.0.8',
            'module load gatk/4.4.0.0',
            'conda activate /share/pkg.8/gatk/4.4.0.0/install/gatk-4.4.0.0',
            cmd,sep = '\n'))
}
