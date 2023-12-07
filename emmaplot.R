


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
      return(p+ggtitle(paste('Enriched pathways for selected genes')))
      
    }else{
      return(p+ggtitle(paste('Enriched pathways with', show_pathway_of)))
      
    }
    
  }else{
    return(p)
    
  }
}


#Annexe functions####

GetPathwaysLinks<-function(res_fgsea,
                   pathway_names=NULL, 
                   col.var="NES",
                   show_pathway_of=NULL,
                   min_edge=0.2,
                   label.size=2.5,
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
    pathways_links<-data.table(get.edgelist(g))
    setnames(pathways_links,new = c('pathway1','pathway2'))
    return(pathways_links)
    
  }else{
    stop('only one pathway gave')
    
  }
}
  
