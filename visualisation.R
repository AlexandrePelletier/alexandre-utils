
#required package
require(pheatmap)
require(stringr)

#DIFFERENTIAL EXPRESSION RESULTS####
#HEATMAPS comparing DEGs 
CompDEGs<-function(res_des,
                   group.by,
                   gene_column='gene',
                   FC_column='log2FoldChange',
                   pval_column='padj',
                   col_range=c(-2.5,2.5),
                   show_rownames=TRUE,
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
    res_des1[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<0.05,'*','')))),.SDcols=pval_column]
    
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
                   cluster_cols = T,
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
                  cluster_cols = T,
                  cellwidth =16,
                  
                  fontsize_number = 10))
  
}




#FGSEA RESULTS####
removeRefKey<-function(term_names)str_remove(str_remove(term_names,'GOBP|GOCC|GOMF|KEGG|NABA|REACTOME|WP|BIOCARTA|PID'),'_')


#EMMAPLOT (pathways network)

{
  require(igraph)
  require(ggraph)
  
#Useful Functions####
#get sim matrix

LeadingEdges<-function(res_fgsea){
  l_genes<-str_extract_all(res_fgsea$leadingEdge,'[A-Za-z0-9]+')
  l_genes<-lapply(l_genes, function(x)x[x!='c'])
  names(l_genes)<-res_fgsea$pathway
  return(l_genes)
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

get_similarity_matrix <- function(leading_edge_list) {
  
  n <- length(leading_edge_list)
  ids<-names(leading_edge_list)
  w <- matrix(NA, nrow=n, ncol=n)
  colnames(w) <- rownames(w) <- names(leading_edge_list)
  
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      w[i,j] <- overlap_ratio(leading_edge_list[ids[i]], leading_edge_list[ids[j]])
    }
  }
  return(w)
}


# get graph of sim
get_igraph <- function(res_fgsea, simmat,leading_edge_list,
                       pathway_names, col.var, min_edge) {
  if(any(duplicated(res_fgsea$pathway)))stop('error: duplicated pathways')
  
  wd <- reshape2::melt(simmat[pathway_names,pathway_names])
  wd <- wd[wd[,1] != wd[,2],]
  # remove NA
  wd <- wd[!is.na(wd[,3]),]
  
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  E(g)$width <- sqrt(wd[, 3] * 5) 
  
  
  
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
  
  res_fgseaf<-res_fgsea[V(g)$name,on='pathway']
  #idx <- unlist(sapply(V(g)$name, function(x) which(x == res_fgseaf$pathway)))
  cnt <- sapply(leading_edge_list, length)

  V(g)$size <- cnt[V(g)$name]

  colVar <- as.numeric(as.vector(res_fgseaf[V(g)$name, on='pathway'][,.SD,.SDcols=col.var][[1]]))

    V(g)$colvar <- colVar

  return(g)
}


#plot the graphs
add_category_nodes <- function(p,col.var,cols=cols) {
  locol=cols[1]
  midcol=ifelse(length(cols==3),cols[2],NULL)
  hicol=cols[-1]
  
  p<-p + ggnewscale::new_scale_fill() +geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ colvar,
                                size =~ size)) +
    scale_size_continuous(name = "number of genes",
                          range = c(3, 8) )
  
  if(!is.null(midcol))p<-p+scale_fill_gradient2(low = locol,mid=midcol, high = hicol,name=col.var,
                                               guide = guide_colorbar()) 
  else p<-p+scale_fill_continuous(low = locol, high = hicol,name=col.var,
                                 guide = guide_colorbar())
  
  p<-p+theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10)) +
    theme(panel.background = element_blank()) 
  return(p)
}
add_node_label <- function(p,label.size=label.size,max.overlaps=10) {
  
  p <- p + geom_node_text(aes_(label=~name),
                          size = label.size, repel=TRUE,
                          max.overlaps=max.overlaps)
  
  return(p)
}

#main function####
emmaplot<-function(res_fgsea,
                   pathway_names=NULL, 
                   col.var="NES",
                   show_pathway_of=NULL,
                   min_edge=0.2,
                   label.size=2.5,
                   cols=c('blue','white','red'),
                   max.overlaps=10){
  require('ggrepel')
  
  if(is.null(pathway_names))pathway_names=res_fgsea[order(pval)]$pathway
  
  lelist<-LeadingEdges(res_fgsea[pathway%in%pathway_names])
  
  if(!is.null(show_pathway_of)){
    
    lelist<-lelist[sapply(lelist, function(leadingedges)any(show_pathway_of%in%leadingedges))]
    if(length(lelist)>0){
      pathway_names<-names(lelist)
      
    }else{
      stop('This gene are not found in any leading edges of the given pathways')
    }
  }
  
  if(length(lelist)>1){
    
    simat<-get_similarity_matrix(lelist)
    
    g <- get_igraph(res_fgsea = res_fgsea,
                    pathway_names = pathway_names,
                    simmat = simat,
                    leading_edge_list = lelist,
                    min_edge = min_edge,
                    col.var = col.var
    )
    
    
    p <- ggraph(g, layout='nicely')
    
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
    ## add dot
    p <- add_category_nodes(p = p,col.var =col.var,cols=cols)
    
    ## add node label
    
    p <- add_node_label(p = p,label.size=label.size,max.overlaps=max.overlaps)
    
  }else{
    p <- ggplot(res_fgsea[pathway%in%pathway_names][,x:=1][,y:=1],aes_string(x='x',y='x'))+
      geom_point(aes_string(size='size',col=col.var))+
      geom_text_repel(aes(label=pathway))+
      scale_color_gradient2(low = cols[1],high = cols[max(1:length(cols))],midpoint = 0,limits=c(-abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var]))),
                                                                                                 abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var])))))+
      theme_graph()
    
  }
 

  if(!is.null(show_pathway_of)){
    if(length(show_pathway_of)>1){
      return(p+ggtitle(paste('Enriched pathways with', show_pathway_of[1],' interactome')))
      
    }else{
      return(p+ggtitle(paste('Enriched pathways with', show_pathway_of)))
      
    }
    
  }else{
    return(p)
    
  }
}





}

#HEATMAPs comparing pathhways###

CompPathways<-function(res_gsea,group.by,legend.compa=NULL,rm.refkey=TRUE,save.pdf=NULL,width =7,height = 7){
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
  res_gsea1[,padjsig:=ifelse(padj<0.001,'***',ifelse(padj<0.01,'**',ifelse(padj<0.05,'*',ifelse(padj<0.25,'.',''))))]

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
  res_de_p<-merge(res_de,degs_pathways,by='gene')
  res_de_p<-unique(res_de_p,by=c('gene','pathway'))
  #create heatmaps
  dep_mat<-data.frame(dcast(res_de_p,gene~pathway,value.var =FC_col),row.names = 'gene')
  dep_mat[is.na(dep_mat)]<-0
  
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
                   # cellwidth =20,
                   # cellheight =  8,
                   
                   fontsize_number = 8))
  }else{
    print(pheatmap(dep_mat,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize= 7,
                   main='Top DEGs',
                   show_rownames = show_rownames,
                   # display_numbers = dep_matp,
                   # cellwidth =20,
                   # cellheight =  8,
                   
                   fontsize_number = 8))
  }
  
  
  
}


