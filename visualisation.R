
#required package
require(pheatmap)
require(stringr)
require(ggplot2)
require(data.table)

source<-function(file,chdir=TRUE)base::source(file,chdir = chdir)

#CompZ: create heatmap with zscore comparing different group for same covariates
CompZ<-function(x,group.by,covcol='cov',zcol='zscore',pvalcol='padj',
                colors=c("blue3", "white", "red3"),
                colorlim=10,color_range=NULL,
                legend.compa=NULL,
                annotation_colors  = NULL,
                cluster_cols=TRUE,
                fontsize_row = 7,cellwidth=12,
                fontsize_number = 12,
                show_rownames=TRUE,
                show_pval=TRUE,
                title=NA,
                FDR_thr=0.05,suggestive_thr=NULL,
                show_several_stars=FALSE,...){
  x<-copy(x)
  
  if('cov'%in%colnames(x)&covcol!='cov')x<-x[,-'cov']
  x[,cov:=.SD,.SDcols=covcol]
  x[,zscore:=.SD,.SDcols=zcol]
  x[,group:=.SD,.SDcols=group.by]
  
  
  x[,zscore:=ifelse(abs(zscore)>colorlim,sign(zscore)*colorlim,zscore)]
  
  matz<-dcast(x,cov~group,value.var ='zscore')|>data.frame(row.names = 'cov')
  matz[is.na(matz)]<-0
  
  if(show_pval){
    if(show_several_stars){
      x[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<=FDR_thr,'*','')))),.SDcols=pvalcol]
    }else{
      x[,padjsig:=lapply(.SD,function(x)ifelse(x<=FDR_thr,'*','')),.SDcols=pvalcol]
    }
    if(!is.null(suggestive_thr)){
      x[,padjsig:=ifelse(.SD[[1]]<=suggestive_thr&.SD[[1]]>FDR_thr,'.',padjsig),.SDcols=pvalcol]
      
    }
 
    matp<-dcast(x,cov~group,value.var ='padjsig')|>data.frame(row.names = 'cov')
    
  }else{
    matp<-FALSE
  }
  

  color_gradient <- colorRampPalette(colors)
  if(is.null(color_range)){
    color_range<-c(-max(abs(matz),na.rm = T),max(abs(matz),na.rm = T))
    
  }else{
    mat_dep<-TRUE
    
  }
  
  if(!is.null(legend.compa)){
    cols_mtd<-c('group',union(group.by,legend.compa))
    mtd_compa<-unique(x[,.SD,.SDcols=cols_mtd])
    mtd_compa[,group:=make.names(group)]
    mtd_compa<-data.frame(mtd_compa,row.names = 'group')[,legend.compa,drop=F]
  }else{
    mtd_compa<-NA
  }
  
  
  col_breaks<-c(((color_range[1]*100):(color_range[2]*100))/100)
  colors=color_gradient(length(col_breaks)-1)
  
  return(pheatmap::pheatmap(matz,
                            breaks =col_breaks,
                            color=colors,
                            fontsize_row = fontsize_row ,
                            show_rownames=show_rownames,
                            main=title,
                            display_numbers = matp,
                            annotation_col = mtd_compa,
                            annotation_colors = annotation_colors,
                            cluster_cols = cluster_cols,
                            number_color = 'white',
                            cellwidth =cellwidth,
                            fontsize_number = fontsize_number,treeheight_row =20,treeheight_col = 20 ,...))
  
  
  
}

CompZ2<-function(x,group.by,covcol='cov',zcol='zscore',pvalcol='padj',
                colors=c("blue3", "white", "red3"),
                colorlim=10,color_range=NULL,
                legend.compa=NULL,
                annotation_colors  = NULL,
                cluster_cols=TRUE,
                fontsize_row = 7,cellwidth=12,
                fontsize_number = 12,
                show_rownames=TRUE,
                show_pval=TRUE,
                title=NULL,
                FDR_thr=0.05,suggestive_thr=NULL,
                show_several_stars=FALSE,...){
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  x<-copy(x)
  
  if(is.null(title)){
    title=zcol
  }
  
  if('cov'%in%colnames(x)&covcol!='cov')x<-x[,-'cov']
  x[,cov:=.SD,.SDcols=covcol]
  x[,zscore:=.SD,.SDcols=zcol]
  x[,group:=.SD,.SDcols=group.by]
  
  
  x[,zscore:=ifelse(abs(zscore)>colorlim,sign(zscore)*colorlim,zscore)]
  
  matz<-dcast(x,cov~group,value.var ='zscore')|>data.frame(row.names = 'cov')
  matz[is.na(matz)]<-0
  
  if(show_pval){
    if(show_several_stars){
      x[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<=FDR_thr,'*','')))),.SDcols=pvalcol]
    }else{
      x[,padjsig:=lapply(.SD,function(x)ifelse(x<=FDR_thr,'*','')),.SDcols=pvalcol]
    }
    if(!is.null(suggestive_thr)){
      x[,padjsig:=ifelse(.SD[[1]]<=suggestive_thr&.SD[[1]]>FDR_thr,'.',padjsig),.SDcols=pvalcol]
      
    }
    
    matp<-dcast(x,cov~group,value.var ='padjsig')|>data.frame(row.names = 'cov')|>as.matrix()
    
  }
  
  
  color_gradient <- colorRampPalette(colors)
  if(is.null(color_range)){
    color_range<-c(-max(abs(matz),na.rm = T),max(abs(matz),na.rm = T))
    
  }else{
    mat_dep<-TRUE
    
  }
  
  if(!is.null(legend.compa)){
    cols_mtd<-c('group',union(group.by,legend.compa))
    mtd_compa<-unique(x[,.SD,.SDcols=cols_mtd])
    mtd_compa[,group:=make.names(group)]
    mtd_compa<-data.frame(mtd_compa,row.names = 'group')[colnames(matz),legend.compa,drop=F]
    if(is.null(annotation_colors)){
      ha <- HeatmapAnnotation(
        df = mtd_compa
      )
    }else{
      ha <- HeatmapAnnotation(
        df = mtd_compa,
        col = annotation_colors
      )
    }

  }else{
    ha<-NULL
  }
  

  
  col_breaks<-c(((color_range[1]*100):(color_range[2]*100))/100)
  colors=color_gradient(length(col_breaks))
  


# ---- Color scale (breaks + colors) ----
col_fun <- circlize::colorRamp2(col_breaks, colors)



# ---- Heatmap ----
ht <- Heatmap(
  matz,
  col = col_fun,
  
  # titles
  name = title,
  
  # clustering
  cluster_columns = cluster_cols,
  
  # row names
  show_row_names = show_rownames,
  row_names_gp = gpar(fontsize = fontsize_row),
  
  # dendrogram size
  row_dend_width = unit(20, "mm"),
  column_dend_height = unit(20, "mm"),
  
  # annotation
  top_annotation = ha,
  
  # cell size (approximate)
  width = ncol(matz) * unit(cellwidth, "points"),
  
  # ---- numbers in cells ----
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- matp[i, j]
    
    if (is.na(val)) return()
    
    txt <- if (is.numeric(val)) {
      sprintf("%.2f", val)
    } else {
      as.character(val)
    }
    
    grid.text(
      txt,
      x, y,
      gp = gpar(
        fontsize = ifelse(str_detect(val,'\\*'),fontsize_number,fontsize_number+8),
        col = "white"
      )
    )
  }
)

return(draw(ht))
}

#DIFFERENTIAL EXPRESSION RESULTS####
#CompDEGs: HEATMAPS comparing DEGs of different groups/comparison
#INPUT:
#res_des: data.table of the differential expression results for the genes to compare. 
#group.by: the column(s) in this data.table of the groups to compare
#gene_column: name of the column with the gene name of the DEGs, 
#FC_column: name of the column with the log2 Fold change of the DEGs, 
#pval_column: name of the column with the (adjuster) pvalue of the DEGs,
#col_range: adjust the minimum and maximum value of the fold change colors
#save.pdf: where to save the heatmap as pdf. default NULL (not save as pdf)
#OUTPUT: a pheatmap heatmap with asterisk for each gene-comparison if p<0.001 : ***, p<0.01: **, p<0.05: *, p<0.25: '.'. 
CompDEGs<-function(res_des,
                   group.by,
                   gene_column='gene',
                   FC_column='log2FoldChange',
                   pval_column='padj',
                   col_range=c(-2.5,2.5),
                   colors=c("blue3", "white", "red3"),
                   colors_resol=100,
                   na_value=0,
                   show_NA=TRUE,
                   show_rownames=TRUE,
                   cluster_cols=TRUE,
                   show_pval=TRUE,
                   asterisk_fontsize = 10,
                   asterisk_color = 'white',
                   width =7,
                   height = 7){
  require('pheatmap')
  require('data.table')
  
  res_des1<-copy(res_des)
  if(length(group.by)>1){
    
  }
  
  res_des1[,comparison:=apply(.SD,1,function(...)paste(...,collapse = '.')),.SDcols=group.by]
  
  res_des1[,gene:=.SD,.SDcols=gene_column]
  
  mat_de<-data.frame(dcast(res_des1,
                             gene~comparison,value.var =FC_column),row.names = 'gene')
  
  mat_de[is.na(mat_de)]<-na_value
  if(show_pval){
    #add pvalue
    res_des1[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<0.05,'*',ifelse(x<0.25,'.',''))))),.SDcols=pval_column]
    
    mat_dep<-data.frame(dcast(res_des1,gene~comparison,value.var ='padjsig'),row.names = 'gene')
    if(!show_NA){
      mat_dep[is.na(mat_dep)]<-''
      
    }
    
  }else{
    mat_dep<-FALSE
  }
  
  color_gradient <- colorRampPalette(colors)
  
  
  col_breaks<-c(((col_range[1]*colors_resol):(col_range[2]*colors_resol))/colors_resol)
  colors=color_gradient(length(col_breaks)-1)

  return(pheatmap::pheatmap(mat_de,
                  breaks =col_breaks,
                  color=colors,
                  fontsize_row = 7,
                  show_rownames=show_rownames,
                  main=FC_column,
                  display_numbers = mat_dep,
                  cluster_cols = cluster_cols,
                  cellwidth =16,
                  number_color = asterisk_color,
                  fontsize_number = asterisk_fontsize))
  
}




#FGSEA RESULTS####
removeRefKey<-function(term_names)str_remove(str_remove(term_names,'GOBP|GOCC|GOMF|KEGG|NABA|REACTOME|WP|BIOCARTA|PID'),'_')


#EMMAPLOT (pathways network)
source('emmaplot.R')





#HEATMAPs comparing pathhways###

#legend.compa /.pathway: colomn name in the dataframe to annotate the variable or the pathways
CompPathways<-function(res_gsea_or_or,group.by,
                       legend.compa=NULL,
                       legend.pathway=NULL,
                       rm.refkey=FALSE,
                       show_pval=TRUE,
                       pval_col='padj',effect_col='NES',pathw_col='pathway',
                       padj_sig_thr=0.1,
                       padj_sugg_thr=0.25,
                       signif_col=NULL,
                       colors=c("blue3", "white", "red3"),
                       white_thr=0.5,
                       colors_resol=100,
                       na_value=0,
                       show_NA=TRUE,
                       width =7,height = 7,max_color=2,
                       cellwidth=16,cluster_cols=TRUE,
                       cluster_rows=TRUE,main=NULL,fontsize_row=7,fontsize_number=10,...){
  require('pheatmap')
  require('data.table')
  res_gsea1<-copy(res_gsea_or_or)
  
  if(is.null(main)){
    main=effect_col
  }
  
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
  
  mat_gsea[is.na(mat_gsea)]<-na_value


  if(show_pval){
    #add pvalue
    res_gsea1[,padj:=.SD,.SDcols=pval_col]
    if(is.null(signif_col)){
      res_gsea1[,padjsig:=ifelse(padj<0.001,'***',ifelse(padj<0.01,'**',ifelse(padj<padj_sig_thr,'*',ifelse(padj<padj_sugg_thr,'.',''))))]
      
    }else{
      res_gsea1[,padjsig:=.SD,.SDcols=signif_col]
      
    }
    res_gsea1[is.na(padj),padjsig:='']
    mat_gseap<-data.frame(dcast(res_gsea1,pathw~comp,value.var ='padjsig'),row.names = 'pathw')[rownames(mat_gsea),colnames(mat_gsea)]
    
    if(!show_NA){
      mat_gseap[is.na(mat_gseap)]<-''
      
    }
  }else{
    mat_gseap<-FALSE
  }
  
  
  color_gradient <- colorRampPalette(colors)
  
  if(all(c(-1,1)%in%unique(sign(res_gsea1[[effect_col]])))){
    #if(is.null(revert_color))revert_color=TRUE
    
    col_breaks<-c((-(colors_resol*max_color):(colors_resol*max_color))/colors_resol)
    col_breaks<-col_breaks[col_breaks>white_thr|col_breaks<(-white_thr)]
    #if(is.null(palette_name))palette_name='RdBu'

  
  }else{
    #if(is.null(revert_color))revert_color=FALSE
    col_breaks<-0:(colors_resol*max_color)/colors_resol
    col_breaks<-col_breaks[col_breaks>white_thr|col_breaks<(-white_thr)]
    
    #if(is.null(palette_name))palette_name='Reds'

  }
  
  # colors=RColorBrewer::brewer.pal(n = 6, name =
  #                                   palette_name)
  
  colors=color_gradient(length(col_breaks)-1)
  # if(revert_color){
  #   colors=rev(colors)
  #   
  # }
  # colors=colorRampPalette(colors)(length(col_breaks)-1)
  
 
  if(!is.null(legend.compa)){
    cols_mtd<-c('comp',union(group.by,legend.compa))
    mtd_compa<-unique(res_gsea1[,.SD,.SDcols=cols_mtd])
    mtd_compa[,comp:=make.names(comp)]
    mtd_compa<-data.frame(mtd_compa,row.names = 'comp')[,legend.compa,drop=F]
  }else{
    mtd_compa<-NA
  }
  
  if(!is.null(legend.pathway)){
    
    cols_mtd<-c('pathw',legend.pathway)
    mtd_path<-unique(res_gsea1[,.SD,.SDcols=cols_mtd],by='pathw')
    # mtd_path[,pathw:=make.names(pathw)]
    mtd_path<-data.frame(mtd_path,row.names = 'pathw')[,legend.pathway,drop=FALSE]

  }else{
    mtd_path<-NA
  }

  return(pheatmap::pheatmap(mat_gsea,
                 breaks =col_breaks,
                 color=colors,
                 fontsize_row = fontsize_row,
                 main=main,
                 display_numbers = mat_gseap,
                 number_color = 'white',
                 cluster_cols = cluster_cols,
                 cellwidth =cellwidth,
                 annotation_col = mtd_compa,
                 annotation_row = mtd_path,
                 fontsize_number = fontsize_number,
                 cluster_rows = cluster_rows,...))
  
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
    print(pheatmap::pheatmap(mat_gsea,
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
  return(pheatmap::pheatmap(mat_gsea,
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
    print(pheatmap::pheatmap(mat_enr,
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
  return(pheatmap::pheatmap(mat_enr,
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
    print(pheatmap::pheatmap(dep_mat,
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
    print(pheatmap::pheatmap(dep_mat,
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

#ModulesHeatmap####
#gene_modules/cell_module: module for  gene i /cell j of mat i*j, will order the modules according to
#module_colors: named vector of colors for each modules
#anno_genes/anno_cells : data.frame of annot for each genes/ cells
#col_genes/col_cells : named list of colors for each annot. (colors have to be named, for categorical) 
ModulesHeatmap<-function(mat,gene_modules,cell_modules,genes_to_label,
                         module_title='Module',module_colors=NULL,
                         anno_cells=NULL,col_cells=NULL,
                         anno_genes=NULL,col_genes=NULL,
                         clamp=c(-2,3),width=5,height=6,fontsize=7){
  require(ComplexHeatmap)
  ## Clamp pour rester dans [-2, 3]
  mat <- pmax(pmin(mat, clamp[2]), clamp[1])
  
  # ── 3. ANNOTATIONS ──────────────────────────────────────────
  n_genes=nrow(mat)
  n_cells=ncol(mat)
  n_modules=length(table(gene_modules))

  all_module_levels <- sort(union(gene_modules, cell_modules))
  
  if(all(str_detect(all_module_levels,'[0-9]'))){
    all_module_levels<-all_module_levels[order(as.numeric(str_extract(all_module_levels,'[0-9]+')))]
    
  }
  
  gene_modules <- factor(gene_modules, levels = all_module_levels)
  cell_modules <- factor(cell_modules, levels = all_module_levels)

  ## Palette modules
  if(is.null(module_colors)){
    module_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(length(all_module_levels), "Paired"))(length(all_module_levels)),
      all_module_levels
    )
  }
  
  

  if(is.null(anno_cells)){
    anno_cells=data.frame(foo=1:length(cell_modules))
  }
  anno_cells[[module_title]]=cell_modules
  

  if('foo'%in%colnames(anno_cells)){
    anno_cells[['foo']]<-NULL
  }
  
  
  if(is.null(col_cells)){
    col_cells=list()
  }
  
  col_cells[[module_title]]=module_colors
  
  ## Annotation colonnes (cellules)
  col_ann <- HeatmapAnnotation(
    df = anno_cells,
    col    = col_cells,
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  
  
  if(is.null(anno_genes)){
    anno_genes=data.frame(foo=1:length(gene_modules))
  }
  anno_genes[[module_title]]=gene_modules
  
  if('foo'%in%colnames(anno_genes)){
    anno_genes[['foo']]<-NULL
  }
  
  if(is.null(col_genes)){
    col_genes=list()
  }
  
  col_genes[[module_title]]=module_colors
  
  ## Annotation lignes (gènes)
  row_ann <- rowAnnotation(
    df = anno_genes,
    col    = col_genes,
    show_legend = FALSE
  )
  
  # ── 4. LABELS DE GÈNES SÉLECTIONNÉS ─────────────────────────
  
  labeled_genes <- unlist(genes_to_label)
  
  ## Position de chaque gène labellé dans la matrice
  label_idx <- match(labeled_genes, rownames(mat))
  label_idx <- label_idx[!is.na(label_idx)]
  
  ## Annotation texte sur la droite
  gene_label_ann <- rowAnnotation(
    Genes = anno_mark(
      at    = label_idx,
      labels = rownames(mat)[label_idx],
      labels_gp = gpar(fontsize = fontsize, fontface = "italic"),
      link_width = unit(4, "mm"),
      padding    = unit(1, "mm")
    )
  )
  
  # ── 5. PALETTE DE COULEURS ───────────────────────────────────
  
  col_fun <- circlize::colorRamp2(
    breaks = c(-2, -0.5, 0, 0.5, 1.5, 3),
    colors = c("#2166AC", "#92C5DE", "#F7F7F7", "#FDDBC7", "#D6604D", "#B2182B")
  )
  
  # ── 6. DESSIN DE LA HEATMAP ─────────────────────────────────
  
  ht <- Heatmap(
    mat,
    name = "Scaled\nexpression",
    
    ## Couleurs
    col = col_fun,
    
    ## Annotations
    top_annotation   = col_ann,
    left_annotation  = row_ann,
    right_annotation = gene_label_ann,
    
    ## Pas de clustering (ordre par module voulu)
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    
    ## Séparateurs entre modules
    row_split    = gene_modules,
    column_split = cell_modules,
    row_gap      = unit(1, "mm"),
    column_gap   = unit(1, "mm"),
    
    ## Titres de blocs
    row_title_gp    = gpar(fontsize = fontsize+1, fontface = "bold"),
    column_title_gp = gpar(fontsize = fontsize+1, fontface = "bold"),
    row_title_rot   = 0,
    
    ## Noms de gènes / cellules masqués (trop nombreux)
    show_row_names    = FALSE,
    show_column_names = FALSE,
    
    ## Légende
    heatmap_legend_param = list(
      title_gp    = gpar(fontsize = fontsize+2, fontface = "bold"),
      labels_gp   = gpar(fontsize = fontsize+2),
      legend_height = unit(0.5, "in")
    ),
    
    ## Taille de la figure
    width  = unit(width, "in"),
    height = unit(height, "in")
  )
  
  # ── 7. EXPORT ────────────────────────────────────────────────
  
  return(draw(ht,
              heatmap_legend_side = "right",
              annotation_legend_side = "right",
              merge_legend = TRUE
  ))
}


