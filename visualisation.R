
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

  
  wd <- reshape2::melt(simmat[pathway_names,pathway_names])
  wd <- wd[wd[,1] != wd[,2],]
  # remove NA
  wd <- wd[!is.na(wd[,3]),]
  
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  E(g)$width <- sqrt(wd[, 3] * 5) 
  
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
  idx <- unlist(sapply(V(g)$name, function(x) which(x == res_fgsea$pathway)))
  cnt <- sapply(leading_edge_list[idx], length)
  
  V(g)$size <- cnt
  
  colVar <- as.numeric(as.vector(res_fgsea[idx, ..col.var][[1]]))
  V(g)$colvar <- colVar
  return(g)
}


#plot the graphs
add_category_nodes <- function(p,col.var) {
  
  p<-p + ggnewscale::new_scale_fill() +
    geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ colvar,
                                size =~ size)) +
    scale_size_continuous(name = "number of genes",
                          range = c(3, 8) ) +
    scale_fill_continuous(low = "blue", high = "red",name=col.var,
                          guide = guide_colorbar()) + 
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10)) +
    theme(panel.background = element_blank()) 
  return(p)
}
add_node_label <- function(p) {
  
  p <- p + geom_node_text(aes_(label=~name),
                          size = 2.5, repel=TRUE)
  
  return(p)
}

#main function####
emmaplot<-function(res_fgsea,
                   pathway_names, 
                   col.var="NES",
                   min_edge=0.2){
  
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
  p <- add_category_nodes(p = p,col.var =col.var)
  ## add node label
  
  p <- add_node_label(p = p)
  
  
  return(p)
}





}


