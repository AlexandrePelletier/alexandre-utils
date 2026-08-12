

library(igraph)
library(ggraph)
library(stringr)
library(ggplot2)
library(data.table)


#Useful Functions####
#get sim matrix

LeadingEdges<-function(res_fgsea){
    if(is.character(res_fgsea$leadingEdge)){
      l_genes<-strsplit(res_fgsea$leadingEdge,'\\|')
    }else{
      l_genes<-res_fgsea$leadingEdge
    }
      

    names(l_genes)<-res_fgsea$pathway
    return(l_genes)
}

  
 


overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}


get_similarity_matrix <- function(genesets) {
  # w=sapply(genesets, function(x)sapply(genesets, function(y)length(intersect(x,y))/length(union(x,y))))
  # w[lower.tri(w,diag = TRUE)]=NA
  # 
  require(Matrix)
  
  # Example:
  # genesets is a named list
  # genesets <- list(
  #   GS1 = c("A","B","C"),
  #   GS2 = c("B","C","D"),
  #   GS3 = c("X","Y")
  # )
  
  # ---- Build sparse incidence matrix ----
  
  all_genes <- unique(unlist(genesets))
  
  i <- rep(seq_along(genesets), lengths(genesets))
  j <- match(unlist(genesets), all_genes)
  
  M <- sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(length(genesets), length(all_genes)),
    dimnames = list(names(genesets), all_genes)
  )
  
  # ---- Pairwise intersections ----
  # tcrossprod gives shared gene counts
  
  intersections <- tcrossprod(M)
  
  # ---- Set sizes ----
  
  sizes <- Matrix::rowSums(M)
  
  # ---- Pairwise unions ----
  
  unions <- outer(sizes, sizes, "+")
  
  # unions = |A| + |B| - |A∩B|
  unions <- unions - as.matrix(intersections)
  
  # ---- Jaccard similarity ----
  
  jaccard <- as.matrix(intersections) / unions
  
  jaccard[lower.tri(jaccard,diag = TRUE)] <-NA
  
  return(jaccard)
}
 

# get graph of sim
get_igraph <- function(res_fgsea, simmat,genesets,node_size=NULL,
                       pathway_names=NULL, min_edge=0.2, col.var=NULL) {
  if(any(duplicated(res_fgsea$pathway)))stop('error: duplicated pathways')
  
  if(is.null(pathway_names)){
    pathway_names<-res_fgsea$pathway
  }
  wd <- reshape2::melt(simmat[pathway_names,pathway_names])
  wd <- wd[wd[,1] != wd[,2],]
  # remove NA
  wd <- wd[!is.na(wd[,3]),]
  
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  
  E(g)$width <- wd[, 3] # * 5*0.2/min_edge
  
  
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- igraph::delete.edges(g, E(g)[wd[, 3] < min_edge])

  res_fgseaf<-res_fgsea[V(g)$name,on='pathway']
  #idx <- unlist(sapply(V(g)$name, function(x) which(x == res_fgseaf$pathway)))
  if(is.null(node_size)){
    cnt <- sapply(genesets, length)
    
    V(g)$size <- cnt[V(g)$name]
  }else{
    V(g)$size <- node_size[V(g)$name]
    
  }

  
  if(!is.null(col.var)){
    colVar <- res_fgseaf[V(g)$name, on='pathway'][[col.var]]
    V(g)$colvar <- colVar
    
  }

  return(g)
}


#plot the graphs
add_category_nodes <- function(p,col.var,cols=c('blue','white','red'),cols_lims=NULL,node_size_name="number of genes") {
  
  if(!any('discrete'%in%cols)){
    locol=cols[1]
    if(length(cols)==3)midcol=cols[2]
    
    hicol=ifelse(length(cols)==3,cols[3],cols[2])
    
  }

  p<-p + ggnewscale::new_scale_fill() +geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ colvar,
                                                                   size =~ size)) +
    
    scale_size_continuous(name = node_size_name,
                          range = c(3, 8) )
  
  if('discrete'%in%cols){
    p<-p+scale_fill_discrete(name=col.var ) 
    
  }else if(length(cols)==3){
    p<-p+scale_fill_gradient2(low = locol,mid=midcol, high = hicol,name=col.var,
                              guide = guide_colorbar(),
                              limits=cols_lims,midpoint = 0) 
    
  }else{
    p<-p+scale_fill_continuous(low = locol, high = hicol,name=col.var,
                               guide = guide_colorbar(),
                             limits=cols_lims)


  }
  
  p<-p+theme(legend.title = element_text(size = 10),
             legend.text  = element_text(size = 10)) +
    theme(panel.background = element_blank()) 
  
  return(p)
}
add_node_label <- function(p,label.size=label.size,max.overlaps=10) {

  p <- p + geom_node_text(aes(label=name),
                          size = label.size, repel=TRUE,
                          max.overlaps=max.overlaps)
  
  return(p)
}


FormatEnrichmentRes<-function(x,pathway_col=NULL){
  x<-copy(x)
  if(!is.null(pathway_col)){
    
    x[['pathway']]<-x[[pathway_col]]
  }
  
  if('n.overlap'%in%colnames(x)){
    if(is.null(pathway_col)){

      x[,pathway:=term]
      
    }
    x[,size:=n.overlap]
    if('genes.overlap'%in%colnames(x)){
      x[,leadingEdge:=genes.overlap]
    }else{
      warning('genes.overlap column not present while ORA format')
    }
    
  }
  return(x)
  
}

#main function####
#EMMAPLOT
#col.var: color nodes based on which variable in the res_fgsea
#size.var: size of the nodes based on which variable in the res_fgsea
#show_pathway_of: vector of gene. Show only pathway enriched containing this gene(s)
#simat: the similarity/correlation matrix
#min_edge: score from the simat from which to create edge

emmaplot<-function(res_fgsea,
                   pathway_names=NULL, 
                   pathway_col=NULL,
                   col.var=NULL,
                   size.var=NULL,
                   group.by=NULL,
                   show_pathway_of=NULL,
                   simat=NULL,
                   min_edge=0.2,
                   label.size=2.5,
                   cols=c('blue3','white','red3'),
                   max.overlaps=10,cols_lims=NULL){
  if(!is.null(group.by)){
    
    ps<-lapply(split(res_fgsea,by=group.by), function(res)emmaplot(res,  pathway_names=pathway_names, 
                                                                   pathway_col=pathway_col,
                                                                   col.var=col.var,
                                                                   size.var=size.var,
                                                                   show_pathway_of=show_pathway_of,
                                                                   simat=simat,
                                                                   min_edge=min_edge,
                                                                   label.size=label.size,
                                                                   cols=cols,
                                                                   max.overlaps=10,cols_lims=cols_lims))
    ps<-lapply(names(ps), function(n)ps[[n]]+ggtitle(n))
    return(wrap_plots(ps))
  }
  require('ggrepel')
  if(!'data.table'%in%class(res_fgsea)){
    res_fgsea<-data.table(res_fgsea)
  }else{
    res_fgsea<-copy(res_fgsea)
    
  }
  if(!is.null(pathway_col)){
    res_fgsea<-copy(res_fgsea[,pathway:=.SD,.SDcols=pathway_col])
  }
  
  if(!is.null(cols_lims)&length(cols)>2){
    
    if(!(any(cols_lims<0)&any(any(cols_lims>=0)))){
      cols=c('white','red')
    }
  }


  res_fgsea<-FormatEnrichmentRes(res_fgsea,pathway_col = pathway_col )
  
  if(is.null(col.var)){
    if(all(c('term','n.overlap')%in%colnames(res_fgsea))){
      col.var='fold.enrichment'
      
    }else{
      col.var='NES'
      
    }
  }

  

 
    if(all(table(res_fgsea[[col.var]])>1)|!is.numeric(res_fgsea[[col.var]])){
      
      res_fgsea[,(col.var):=lapply(.SD,as.character),.SDcols=col.var]
      
      cols='discrete'
    }
  
  
  
  if(is.null(pathway_names))pathway_names=res_fgsea$pathway
  
 
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
    
    if(is.null(simat)){
      simat<-get_similarity_matrix(lelist)
      
    }else{
      simat<-simat[pathway_names,pathway_names]
      diag(simat)<-NA
    }
    
    
    if(!is.null(size.var)){
      node_size=setNames(res_fgsea[[size.var]],res_fgsea$pathway)
      node_size_name=size.var
    }else{
      node_size=NULL
      node_size_name=NULL
    }
    g <- get_igraph(res_fgsea = res_fgsea,
                    pathway_names = pathway_names,
                    simmat = simat,
                    node_size =node_size ,
                    genesets = lelist,
                    min_edge = min_edge,
                    col.var = col.var
    )
    
    
    p <- ggraph(g, layout='nicely')
    #width=enquo(width)
    p <- p + geom_edge_link(alpha=.7, aes(edge_width=width),
                            colour='gray')+scale_edge_width_continuous(range = c(0.5,3),
                                                                           limits = c(min_edge,max(E(g)$width)))
    ## add dot
    p <- add_category_nodes(p = p,col.var =col.var,cols=cols,cols_lims=cols_lims,
                            node_size_name =node_size_name)
    
    ## add node label
    
    p <- add_node_label(p = p,label.size=label.size,max.overlaps=max.overlaps)
    
  }else{
    p <- ggplot(res_fgsea[pathway%in%pathway_names][,x:=1][,y:=1],aes_string(x='x',y='x'))+
      geom_point(aes_string(size=size.var,col=col.var))+
      geom_text_repel(aes(label=pathway))+
      scale_color_gradient2(low = cols[1],high = cols[length(cols)],
                            midpoint = 0,limits=c(-abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var]))),
                                                                                                 abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var])))))+
      theme_void()
    
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
                   show_pathway_of=NULL,
                   min_edge=0.2,pathway_col=NULL){

  res_fgsea<-FormatEnrichmentRes(res_fgsea,pathway_col = pathway_col)
  
  if(is.null(pathway_names))pathway_names=unique(res_fgsea[order(pval)]$pathway)
  
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
                    genesets = lelist,
                    min_edge = min_edge
    )
    pathways_links<-data.table(get.edgelist(g))
    setnames(pathways_links,new = c('pathway1','pathway2'))
    pathways_links[,weight:=simat[pathway1,pathway2],by=c('pathway1','pathway2')]
    
    return(pathways_links[!is.na(weight)][])
    
  }else{
    warning('only one pathway gave')
    return(NULL)
    
  }
}


#ClusterPathways: return the enrichment data.tabl with cluster of pathways added to if
ClusterPathways<-function(x,resolution=1,method='leiden',weights=NULL,min_edge=0.2,pathway_col=NULL){
  require(igraph)
  x<-FormatEnrichmentRes(x,pathway_col=pathway_col)
  pathway_col='pathway'
  if(length(unique(x$pathway))>1){
    if(is.null(weights)){
      if(!'weight'%in%colnames(x)){
        # if(is.null(min_edge)){
        #   min_edge=0.2
        # }
        message('getting pathway links weight')
        links<-GetPathwaysLinks(x,min_edge = min_edge,pathway_col=pathway_col)
        
      }
      
    }
    
    
    graph<-graph.data.frame(data.frame(links),directed = F)
    if(is.null(weights)){
      weights=igraph::edge_attr(graph,'weight')
      
    }
    
    
    if(method=='louvain'){
      cl<-cluster_louvain(graph, weights = weights, resolution = resolution)
      
    }else if (method=='leiden'){
      cl<-cluster_leiden(graph,objective_function = 'modularity', weights = weights, 
                         resolution = resolution,
                         n_iterations = 5)
      
    }else{
      stop('only leiden and louvain method implemented yet')
    }
    
    return(merge(x,data.table(pathway=cl$names,cluster=cl$membership),
                by='pathway',all.x=TRUE)[x$pathway,on='pathway'])
    
  }else{
    warning('only one pathway provided, return NA')
    return(x[,cluster:=NA])
  }

}


#PathwaysSelection: cluster the pathways using emmaplot, pick the first by cluster based on pval, 
#and let add or edit the pathways selected by the user. Return the enrichment data.table with 
#res_enr:
#interactive: add or edit the pathways selected (with +/-)
#add_plot: which additional metric to color by the pathways on the graph. if a plot, just plot next to the graph
#split.by: from which variable to split by the emmaplot to select step by step. NA/0 no spliting
PathwaysSelection<-function(res_enr,interactive=TRUE,add_plot=NULL,plots_widths=c(1,1),auto_select_top_cluster=FALSE,
                            resolution =1,min_edge=0.2,pathway_col=NULL,split.by='cluster',group.by=NULL){
  res_enr<-FormatEnrichmentRes(copy(res_enr),pathway_col = pathway_col)
  
  selection<-c()
  
  if(!'selected'%in%colnames(res_enr)){
    res_enr[,selected:=FALSE]
  }
  
  
  if(is.null(group.by)){
    res_enr[,group:=0]
  }else {
    res_enr[,group:=apply(.SD,1,paste,collapse='.'),.SDcols =group.by ]
    
  }
  
  selection<-res_enr[(selected)]$pathway
  
  
  for(g in unique(res_enr$group)){
    message(g)
    res_enrf<-copy(res_enr[group==g])
    res_enrf<-res_enrf[order(pval)]
    res_enrf[,pathway_num:=1:.N]
    res_enrf[,pathway_with_num:=paste(pathway_num,pathway,sep='-')]
    res_enrf<-ClusterPathways(res_enrf,resolution = resolution ,min_edge = min_edge,pathway_col = 'pathway_with_num')
    res_enrf[,pathway:=str_remove(pathway,'^[0-9]+-')]
    res_enrf[is.na(cluster),cluster:=0]
    
    if(auto_select_top_cluster){
      res_enrf[,selected:=rank(pval)==1,by='cluster']
      selection<-union(selection,res_enrf[(selected)]$pathway)
    }
    # cluster pathways
    if(is.null(split.by)){
      res_enrf[,split:=0]
      
    }else if(split.by=='cluster'){

      res_enrf[,split:=cluster]
      
    }else if(all(split.by%in%colnames(res_enrf))){
      res_enrf[,split:=apply(.SD,1,paste,collapse='.'),.SDcols =split.by ]
      
    }else{
      res_enrf[,split:=split.by]
    }
    
    # #select the top per cluster and all 0 cluster
    # res_enrf[,selected:=rank(pval)==1,by='split']
    #   
    #if interactive, show and let choose
      
      for(cl in unique(res_enrf$split)){
        message('cluster ',cl)
        p1<-emmaplot(res_enrf[split==cl],col.var = 'selected',pathway_col = 'pathway_with_num')
        
        if(!is.null(add_plot)){
          if(is.character(add_plot)){
            add_plot<-emmaplot(res_enrf[split==cl],col.var = add_plot,pathway_col = 'pathway_with_num')
            
          }
          pall<-p1+add_plot
          print(pall+plot_layout(nrow = 1,widths = plots_widths))
          
        }else{
          print(p1)
        }
        if(interactive){
          edit<-readline(prompt = 'pathways selection:')
          
          torm<-str_extract_all(edit,'\\-[0-9]+')|>unlist()|>str_extract('[0-9]+')|>as.numeric()
          toadd<-strsplit(edit,'\\+|\\-')|>unlist()|>as.numeric()
          toadd<-setdiff(toadd,torm)
          
          
          selection<-union(setdiff(selection,res_enrf[pathway_num%in%torm]$pathway),
                           res_enrf[pathway_num%in%toadd]$pathway)
          print(selection)
          
          res_enr[pathway%in%selection,selected:=TRUE]
          res_enrf[pathway%in%selection,selected:=TRUE]
        }

        
      }
    
  }
    
   

  return(selection)
  
}

#return the enrichment data table filitred with your saved and new selected pathways
EditPathwaysSelection<-function(res_enr,saved,group.by=NULL,force_edit=FALSE,filter_only=FALSE){
  
  
  selection.old<-saved$term|>unique()
  if(!force_edit){
    res_enrnew<-res_enr[!query%in%unique(saved$query)]
    
  }else{
    res_enrnew<-res_enr
  }
  
  
  if(nrow(res_enrnew)>0&!filter_only){
    
    res_enrnew[,selected:=term%in%selection.old]
    selection<-PathwaysSelection(res_enrnew,group.by =group.by,split.by = NULL)
    selection<-union(selection.old,selection)
    
    
  }else{
    selection<-selection.old
  }
  
  res_enrf<-res_enr[term%in%selection]
  
  if(length(setdiff(selection,selection.old))>0|force_edit){
    res_enrf[,term_nice:=TidyPathwayNames(term,remove_gs_source = TRUE)]
    hm<-CompPathways(res_enrf,
                     group.by = 'query',pathw_col = 'term_nice',effect_col = 'log2FE',max_color = 3)
    res_enrf[,pval10:=ifelse(pval<1e-10,1e-10,pval)]
    res_enrf<-res_enrf[data.table(term_nice=hm$tree_row$labels[hm$tree_row$order],
                                  num=1:length(hm$tree_row$labels)),on=c('term_nice')]
    res_enrf[,term_num:=paste(num,term_nice)]
    
    
    #further filter
    keepediting=TRUE
    while(keepediting){
      p<-ggplot(res_enrf[padj<0.25])+
        geom_point(aes(x=query.,y=term_num,size=precision*100,col=-log10(pval10)))+
        scale_x_discrete(guide = guide_axis(angle = 60))+theme_bw()+
        scale_color_gradient2(high = 'red3',low = 'blue2')+scale_y_discrete(limits=res_enrf[order(num)]$term_num|>unique())+
        labs(size='% of the gene-set',y='MSigDB GO/CP pathways',x='cell state signature')
      print(p)
      to_rm<-readline(prompt = "to remove (the num, sep by any sep. to stop just leave empty):")|>str_extract_all('[0-9]+')|>unlist()|>as.numeric()
      if(length(to_rm[!is.na(to_rm)])>0){
        res_enrf<-res_enrf[!num%in%to_rm]
        
      }else{
        keepediting=FALSE
        
        res_enrf<-res_enrf[,-c('term_num','num')]
      }
      
    }
    
    
    
  }
  return(res_enrf)
  
}


GetRepresentativePathways<-function(res_fgsea,group.by=NULL,pathway_names=NULL,return_full=TRUE,padj.thr=0.05,resolution=1){
  
  res<-copy(res_fgsea)
  if(max(res$padj)>padj.thr){
    message('removing pathways not passing padj threshold ',padj.thr)
    message(nrow(res[padj<padj.thr]),'/',nrow(res),' pathways conserved')
    res<-res[padj<padj.thr]
    
  }
  
  if(!is.null(pathway_names))
    res<-res[pathway%in%pathway_names]
  
  if(is.null(group.by)){
    group.by='group'
    res<-res[,(group.by):=1]
  }
  
  if(length(group.by)>1){
    res<-res[,group:=apply(.SD,1,function(x)paste(x,collapse = '_')),.SDcols=group.by]
    group.by='group'
  }
  
  res<-rbindlist(lapply(unlist(unique(res[,..group.by])), function(g){
    message('finding representative pathways for ',g)
    res1<-unique(res[g,on=group.by][order(padj)],by='pathway')
    if(nrow(res1)>1){
      pathways_links<-GetPathwaysLinks(res1)
      
      pathways_cluster<-ClusterPathways(pathways_links,
                                        method='louvain',resolution = resolution,weights = NA)
      
      
      pathway_stat<-data.table(pathway=union(pathways_links$pathway1,pathways_links$pathway2))
      
      pathway_stat[,n.link:=nrow(pathways_links[pathway1==pathway|pathway2==pathway]),by='pathway']
      
      pathways_links<-merge(pathways_links,copy(pathways_cluster)[,pathway1:=pathway][,cluster1:=cluster][,.(pathway1,cluster1)])
      pathways_links<-merge(pathways_links,copy(pathways_cluster)[,pathway2:=pathway][,cluster2:=cluster][,.(pathway2,cluster2)],by='pathway2')
      pathway_stat<-merge(pathway_stat,pathways_cluster)
      pathway_stat[,n.link.in.cluster:=nrow(pathways_links[(pathway1==pathway&cluster2==cluster)|(pathway2==pathway&cluster1==cluster)]),by='pathway']
      
      res1<-merge(res1,pathway_stat,all.x = TRUE)
      
      res1[is.na(cluster),cluster:=0]
      
      #central one + the one with bigger NES, padj by cluster
      res1[,top.central:=rank(n.link.in.cluster)==max(rank(n.link.in.cluster)),by='cluster']
      res1[,top.NES:=rank(abs(NES))==max(rank(abs(NES))),by='cluster']
      res1[,top.pval:=rank(pval)==min(rank(pval)),by='cluster']
      
      res1[,tops.cluster:=top.central|top.NES|top.pval,by='cluster']
      return(res1)
      
    }else{
      return(res1)
    }
    
  }))
  
  if(return_full)
    return(res)
  else
    return(res[(tops.cluster)])
  
}


