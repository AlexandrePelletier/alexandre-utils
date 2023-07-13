
#visualisation package
library(pheatmap)

#FGSEA RESULTS####
removeRefKey<-function(term_names)str_remove(str_remove(term_names,'GOBP|GOCC|GOMF|KEGG|NABA|REACTOME|WP|BIOCARTA|PID'),'_')

#EMMAPLOT####

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
                   min_edge=0.2,
                   label.size=2.5,
                   cols=c('blue','white','red'),
                   max.overlaps=10){
  if(is.null(pathway_names))pathway_names=res_fgsea[order(pval)]$pathway
  
  lelist<-LeadingEdges(res_fgsea)
  
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
  
  
  return(p)
}





}




#HEATMAPs comparing pathhways

CompPathways<-function(res_gsea,group.by,rm.refkey=TRUE,save.pdf=NULL,width =7,height = 7){
  require('pheatmap')
  require('data.table')
  
  res_gsea1<-copy(res_gsea)
  
  res_gsea1[,comparison:=.SD,.SDcols=group.by]
  
  if(rm.refkey)
    res_gsea1[,pathw:=removeRefKey(pathway)]
  else
    res_gsea1[,pathw:=pathway]
  
    mat_gsea<-data.frame(dcast(res_gsea1,
                             pathw~comparison,value.var ='NES'),row.names = 'pathw')
  
  
  #add pvalue
  res_gsea1[,padjsig:=ifelse(padj<0.001,'***',ifelse(padj<0.01,'**',ifelse(padj<0.05,'*','')))]

  mat_gseap<-data.frame(dcast(res_gsea1,pathw~comparison,value.var ='padjsig'),row.names = 'pathw')
  
  col_breaks<-c((-30:30)/10)
  col_breaks<-col_breaks[col_breaks>0.5|col_breaks<(-0.5)]
  
  if(!is.null(save.pdf)){
    pdf(fp(out,save.pdf),width =7,height = 7)
    print(pheatmap(mat_gsea,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize_row = 7,
                   main='NES',
                   display_numbers = mat_gseap[rownames(mat_gsea),colnames(mat_gsea)],
                   cluster_cols = T,
                   cellwidth =16,
                   
                   fontsize_number = 10))
    
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
                 
                 fontsize_number = 10))
  
}
