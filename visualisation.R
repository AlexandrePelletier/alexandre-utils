
#required package
require(pheatmap)
require(stringr)
require(ggplot2)
require(data.table)

source<-function(file,chdir=TRUE)base::source(file,chdir = chdir)

#DIFFERENTIAL EXPRESSION RESULTS####
#HEATMAPS comparing DEGs 
CompDEGs<-function(res_des,
                   group.by,
                   gene_column='gene',
                   FC_column='log2FoldChange',
                   pval_column='padj',
                   col_range=c(-2.5,2.5),
                   show_rownames=TRUE,
                   cluster_cols=TRUE,
                   show_pval=TRUE,
                   save.pdf=NULL,
                   width =7,
                   height = 7){
  require('pheatmap')
  require('data.table')
  
  res_des1<-copy(res_des)
  
  res_des1[,comparison:=.SD,.SDcols=group.by]
  
  res_des1[,gene:=.SD,.SDcols=gene_column]
  
  mat_de<-data.frame(dcast(res_des1,
                             gene~comparison,value.var =FC_column),row.names = 'gene')
  
  mat_de[is.na(mat_de)]<-0
  if(show_pval){
    #add pvalue
    res_des1[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<0.05,'*',ifelse(x<0.25,'.',''))))),.SDcols=pval_column]
    
    mat_dep<-data.frame(dcast(res_des1,gene~comparison,value.var ='padjsig'),row.names = 'gene')
    
  }else{
    mat_dep<-FALSE
  }
  
  col_breaks<-c(((col_range[1]*10):(col_range[2]*10))/10)

  if(!is.null(save.pdf)){
    
    pdf(save.pdf,width =width,height = height)
    print(pheatmap(mat_de,
                   breaks =col_breaks,
                   show_rownames=show_rownames,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize_row = 7,
                   main=FC_column,
                   display_numbers = mat_dep,
                   cluster_cols = cluster_cols,
                   cellwidth =16,
                   
                   fontsize_number = 8))
    
    dev.off()
  }
  return(pheatmap(mat_de,
                  breaks =col_breaks,
                  color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                        "RdBu")))(length(col_breaks)-1),
                  fontsize_row = 7,
                  show_rownames=show_rownames,
                  main=FC_column,
                  display_numbers = mat_dep,
                  cluster_cols = cluster_cols,
                  cellwidth =16,
                  
                  fontsize_number = 10))
  
}




#FGSEA RESULTS####
removeRefKey<-function(term_names)str_remove(str_remove(term_names,'GOBP|GOCC|GOMF|KEGG|NABA|REACTOME|WP|BIOCARTA|PID'),'_')


#EMMAPLOT (pathways network)
source('emmaplot.R')





#HEATMAPs comparing pathhways###

CompPathways<-function(res_gsea_or_or,group.by,legend.compa=NULL,rm.refkey=FALSE,
                       pval_col='padj',effect_col='NES',pathw_col='pathway',
                       save.pdf=NULL,width =7,height = 7,max_color=2){
  require('pheatmap')
  require('data.table')
  res_gsea1<-copy(res_gsea_or_or)
  
  if(length(group.by)>1){
    res_gsea1[,comp:=Reduce(function(...)paste(...,sep='_'),.SD),.SDcols=group.by]
    
  }else{
    res_gsea1[,comp:=.SD,.SDcols=group.by]
    
  }
  
  if(rm.refkey)
    res_gsea1[,pathw:=removeRefKey(.SD),.SDcols=pathw_col]
  else
    res_gsea1[,pathw:=.SD,.SDcols=pathw_col]
  
  mat_gsea<-data.frame(dcast(res_gsea1,
                           pathw~comp,value.var =effect_col),row.names = 'pathw')

  
  #add pvalue
  res_gsea1[,padj:=.SD,.SDcols=pval_col]
  res_gsea1[,padjsig:=ifelse(padj<0.001,'***',ifelse(padj<0.01,'**',ifelse(padj<0.05,'*',ifelse(padj<0.25,'.',''))))]

  mat_gseap<-data.frame(dcast(res_gsea1,pathw~comp,value.var ='padjsig'),row.names = 'pathw')
  
  col_breaks<-c((-(10*max_color):(10*max_color))/10)
  col_breaks<-col_breaks[col_breaks>0.5|col_breaks<(-0.5)]
  
  if(!is.null(legend.compa)){
    cols_mtd<-c('comp',union(group.by,legend.compa))
    mtd_compa<-unique(res_gsea1[,.SD,.SDcols=cols_mtd])
    mtd_compa[,comp:=make.names(comp)]
    mtd_compa<-data.frame(mtd_compa,row.names = 'comp')[,legend.compa,drop=F]
  }else{
    mtd_compa<-NA
  }
  
  if(!is.null(save.pdf)){
    pdf(save.pdf,width =width,height = width)
    pheatmap(mat_gsea,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize_row = 7,
                   main=effect_col,
                   display_numbers = mat_gseap[rownames(mat_gsea),colnames(mat_gsea)],
                   cluster_cols = T,
                   cellwidth =16,
                   annotation_col =mtd_compa,
                   
                   
                   fontsize_number = 10)
    
    dev.off()

  }
  return(pheatmap(mat_gsea,
                 breaks =col_breaks,
                 color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name =
                                                                       "RdBu")))(length(col_breaks)-1),
                 fontsize_row = 7,
                 main=effect_col,
                 display_numbers = mat_gseap[rownames(mat_gsea),colnames(mat_gsea)],
                 cluster_cols = T,
                 cellwidth =16,
                 annotation_col = mtd_compa,
                 
                 
                 
                 fontsize_number = 10))
  
}


CompGSEA<-function(res_gsea,group.by,legend.compa=NULL,rm.refkey=TRUE,save.pdf=NULL,width =7,height = 7){
  require('pheatmap')
  require('data.table')
  
  res_gsea1<-copy(res_gsea)
  if(length(group.by)>1){
    res_gsea1[,comp:=Reduce(function(...)paste(...,sep='_'),.SD),.SDcols=group.by]
    
  }else{
    res_gsea1[,comp:=.SD,.SDcols=group.by]
    
  }
  
  if(rm.refkey)
    res_gsea1[,pathw:=removeRefKey(pathway)]
  else
    res_gsea1[,pathw:=pathway]
  
  mat_gsea<-data.frame(dcast(res_gsea1,
                             pathw~comp,value.var ='NES'),row.names = 'pathw')
  
  
  #add pvalue
  res_gsea1[,padjsig:=ifelse(padj<0.001,'***',ifelse(padj<0.01,'**',ifelse(padj<0.05,'*','')))]
  
  mat_gseap<-data.frame(dcast(res_gsea1,pathw~comp,value.var ='padjsig'),row.names = 'pathw')
  
  col_breaks<-c((-30:30)/10)
  col_breaks<-col_breaks[col_breaks>0.5|col_breaks<(-0.5)]
  
  if(!is.null(legend.compa)){
    cols_mtd<-c('comp',union(group.by,legend.compa))
    mtd_compa<-unique(res_gsea1[,.SD,.SDcols=cols_mtd])
    mtd_compa[,comp:=make.names(comp)]
    mtd_compa<-data.frame(mtd_compa,row.names = 'comp')[,legend.compa,drop=F]
  }else{
    mtd_compa<-NA
  }
  
  if(!is.null(save.pdf)){
    pdf(save.pdf,width =7,height = 7)
    print(pheatmap(mat_gsea,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize_row = 7,
                   main='NES',
                   display_numbers = mat_gseap[rownames(mat_gsea),colnames(mat_gsea)],
                   cluster_cols = T,
                   cellwidth =16,
                   annotation_col =mtd_compa,
                   
                   
                   fontsize_number = 10))
    
    dev.off()
    dev.off()
  }
  return(pheatmap(mat_gsea,
                  breaks =col_breaks,
                  color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                        "RdBu")))(length(col_breaks)-1),
                  fontsize_row = 7,
                  main='NES',
                  display_numbers = mat_gseap[rownames(mat_gsea),colnames(mat_gsea)],
                  cluster_cols = T,
                  cellwidth =16,
                  annotation_col = mtd_compa,
                  
                  
                  
                  fontsize_number = 10))
  
}

CompGost<-function(res_enr,group.by,score='precision',col_max=1,
                   legend.compa=NULL,save.pdf=NULL,width =7,height = 7){
  require('pheatmap')
  require('data.table')
  
  res_enr1<-copy(res_enr)
  if(length(group.by)>1){
    res_enr1[,comp:=Reduce(function(...)paste(...,sep='_'),.SD),.SDcols=group.by]
    
  }else{
    res_enr1[,comp:=.SD,.SDcols=group.by]
    
  }
  
  
  res_enr1[,pathw:=term_name]
  
  mat_enr<-data.frame(dcast(res_enr1,
                             pathw~comp,value.var =score),row.names = 'pathw')
  
  
  #add pvalue
  res_enr1[,padjsig:=ifelse(p_value<0.001,'***',ifelse(p_value<0.01,'**',ifelse(p_value<0.05,'*','')))]
  
  mat_enrp<-data.frame(dcast(res_enr1,pathw~comp,value.var ='padjsig'),row.names = 'pathw')
  
  if(any(as.matrix(mat_enr)<0)){
    col_breaks<-unique(c(c(1:8/100,1:15/40,0.4,0.45,0.5)*-col_max,
                         c(1:8/100,1:15/40,0.4,0.45,0.5)*col_max))
  }else{
    
    col_breaks<-unique(c(1:8/100,1:15/40,0.4,0.45,0.5))*col_max
    
  }
  col_breaks<-unique(c(1:8/100,1:15/40,0.4,0.45,0.5))*col_max

  if(!is.null(legend.compa)){
    cols_mtd<-unique(c('comp',union(group.by,legend.compa)))
    mtd_compa<-unique(res_enr1[,.SD,.SDcols=cols_mtd])
    mtd_compa[,comp:=make.names(comp)]
    mtd_compa<-data.frame(mtd_compa,row.names = 'comp')[,legend.compa,drop=F]
  }else{
    mtd_compa<-NA
  }
  
  if(!is.null(save.pdf)){
    pdf(save.pdf,width =7,height = 7)
    print(pheatmap(mat_enr,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdYlBu")))(length(col_breaks)-1),
                   fontsize_row = 7,
                   main=score,
                   display_numbers = mat_enrp[rownames(mat_enr),colnames(mat_enr)],
                   cluster_cols = T,
                   cellwidth =16,
                   annotation_col =mtd_compa,
                   
                   
                   fontsize_number = 10))
    
    dev.off()
    dev.off()
  }
  return(pheatmap(mat_enr,
                  breaks =col_breaks,
                  color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                        "RdYlBu")))(length(col_breaks)-1),
                  fontsize_row = 7,
                  main=score,
                  display_numbers = mat_enrp[rownames(mat_enr),colnames(mat_enr)],
                  cluster_cols = T,
                  cellwidth =16,
                  annotation_col = mtd_compa,
                  
                  
                  
                  fontsize_number = 10))
  
}


#HEATMAPS comparing DEGs of pathways
CompDEGsPathways<-function(res_gsea,
                           res_de,
                           top.n=NULL,
                           gene_col='gene_name',
                           FC_col='log2FoldChange',
                           pval_col='padj',
                           col_range=c(-2.5,2.5),
                           transpose=FALSE,
                           show_rownames=FALSE,
                           show_pval=TRUE){
  
  
  #get leading edges
  degs_list<-LeadingEdges(res_gsea)
  
  #filter degs
  if(!is.null(top.n))
    degs_list<-lapply(degs_list, function(x)head(x,top.n))
  
  #trans in dataframe
  degs_pathways<-Reduce(rbind,lapply(names(degs_list),function(p)data.table(pathway=p,
                                                                            gene=degs_list[[p]])))
  
  
  #merge pathway by degs
  res_de[,gene:=.SD,.SDcols = gene_col]
  res_de_p<-merge(res_de,degs_pathways,by='gene')
  res_de_p<-unique(res_de_p,by=c('gene','pathway'))
  #create heatmaps
  dep_mat<-data.frame(dcast(res_de_p,gene~pathway,value.var =FC_col),row.names = 'gene')
  #dep_mat[is.na(dep_mat)]<-0
  
  #add pvalue
  if(show_pval){
    
    res_de_p[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<0.05,'*','')))),.SDcols=pval_col]
    dep_matp<-data.frame(dcast(res_de_p,gene~pathway,value.var ='padjsig'),row.names = 'gene')
    dep_matp[is.na(dep_matp)]<-''
  }
  
  
  #plot heatmap
  if (transpose) {
    if(show_pval) dep_matp<-t(dep_matp)
    dep_mat<-t(dep_mat)
    
  }
  col_breaks<-c(((col_range[1]*10):(col_range[2]*10))/10)
  if(show_pval){
    print(pheatmap(dep_mat,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize= 7,
                   main='Top DEGs',
                   show_rownames = show_rownames,
                   display_numbers = dep_matp,
                   cluster_rows = F,cluster_cols = F,na_col = 'black',
                   # cellwidth =20,
                   # cellheight =  8,
                   
                   fontsize_number = 8))
  }else{
    print(pheatmap(dep_mat,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize= 7,
                   cluster_rows = F,cluster_cols = F,
                   
                   main='Top DEGs',na_col = 'black',
                   show_rownames = show_rownames,
                   # display_numbers = dep_matp,
                   # cellwidth =20,
                   # cellheight =  8,
                   
                   fontsize_number = 8))
  }
  
  
  
}


