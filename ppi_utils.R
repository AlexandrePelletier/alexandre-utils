#Functions for PPI analysis (StringDB)
#StringDB analysis ####
#ppiNeighbors
#input : genes: vector of genes to get neighboorood genes based on STRINGdb
#combined_score.thr: minimum combined_score.thr to be considered as neighbour 
ppiNeighbors<-function(genes,ppi='/projectnb/tcwlab/RefData/STRINGdb/9606.protein.links.detailed.v12.0_reformat.annot.csv.gz',
                       combined_score.thr=250){
  if(is.character(ppi)){
    ppi<-fread(ppi)
    
  }
  neighbs<-rbindlist(lapply(genes, function(g){
    neighb<-ppi[preferred_name.1==g&combined_score>combined_score.thr]
    neighb<-neighb[,.(protein1,preferred_name.1,protein2,preferred_name.2,combined_score)]
    setnames(neighb,
             c('preferred_name.1','preferred_name.2'),
             c('gene_name','neighbor_name'))
    setnames(neighb,
             c('protein1','protein2'),
             c('protein_id','neighbor_id'))
    return(neighb)
    
  }))
  return(neighbs)
  
}

#PlotString: connect genes based on StringDB, by default if combined_score > 250
#genes: genes to plot network of interaction
#res_de: differential expression results to plot statistic
#
PlotString<-function(genes,ppi='/projectnb/tcwlab/RefData/STRINGdb/9606.protein.links.detailed.v12.0_reformat.annot.csv.gz',
                     combined_score.thr=250,res_de=NULL,
                     stat_col='statistic',gene_col='auto',
                     color_limits=NULL,size_limits=NULL,
                     genes_to_display=NULL){
  require(network)
  require(sna)
  require(ggraph)
  require(ggnetwork)
  
  if(is.character(ppi)){
    ppi<-fread(ppi)
    
  }
  
  if(is.null(genes_to_display)){
    genes_to_display=genes
  }
  
  ne=ppiNeighbors(genes,ppi = ppi,combined_score.thr=combined_score.thr)
  
  
  #filter to get only common genes interaction of the LE (between them)
  nef<-ne[neighbor_name%in%genes]
  
  #rm multiple link
  nef[,link:=unique(paste(sort(c(gene_name,neighbor_name)),collapse = '_')),by=c('neighbor_name','gene_name')]
  
  nef<-nef[!duplicated(link)]
  
  nef<-unique(nef[,.(gene_name,neighbor_name,combined_score)])
  #create networks
  net<-as.network(nef[,.(gene_name,neighbor_name,combined_score)],loops = F,directed = F)
  
  
  #add node annot
  if(!is.null(res_de)){
    if(gene_col=='auto'){
      gene_col=names(which(sapply(res_de,
                                  function(x)sum(network.vertex.names(net)%in%x)>length(network.vertex.names(net))*0.5)))
      if(length(gene_col)!=1){
        stop('cannot find gene column of the differential results automatically, need to specify gene_col')
      }
    }
    net %v% stat_col = res_depf[network.vertex.names(net),on=gene_col][[stat_col]]
    
  }else{
    stat_col=NULL
  }
  
  
  #using ggnetwork
  p<-ggplot(ggnetwork(net,weights = "combined_score"), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size = combined_score), color = "grey",alpha=0.5) +
    geom_nodes(aes_string(color = stat_col),size=7) +
    theme_blank()+
    geom_nodetext_repel(aes(label=ifelse(vertex.names%in%genes_to_display,
                                         vertex.names,'')))+
    scale_size(range = c(0,2.5), limits = size_limits)
  
  if(is.numeric(res_depf[[stat_col]])){
    p<-p+scale_color_gradient2(low = 'blue3',mid = 'white',
                              high='red3',midpoint = 0,
                              limits=color_limits)
    
    
  }else{
    p<-p+scale_color_manual(values = color_limits)
    
    
    
  }
  
  return(p)
}


#PermutString
#find if a connexions between a gene and a group of genes is significant based on permutation

PermutString<-function(gene,modules,ppi,connex_genes=NULL,combined_score.thr=250,n.perm=1000,pct.gap.nconnex=0.1,addNeigbConnex=FALSE,nDTthreads=NULL){
  if(is.null(nDTthreads)){
    setDTthreads(threads = parallel::detectCores()-2)
    
  }else{
    setDTthreads(threads = nDTthreads)
    
  }
  message(getDTthreads(),' threads used for data.table')
  
  
  if(is.null(connex_genes)){
    if('combined_score'%in%colnames(ppi)){
      ppif<-unique(ppi[combined_score>combined_score.thr][,.(preferred_name.1,preferred_name.2)])
      
    }else{
      ppif<-ppi
      
    }
    
    #get genes link to the gene
    
    connex_genes<-split(ppif$preferred_name.2,ppif$preferred_name.1)
    
    if(addNeigbConnex){
      #get the genes connected to the gene of int and its neighbor
      connex_genes<-lapply(connex_genes, function(con_genes){
        return(unique(c(con_genes,ppif[preferred_name.1%in%con_genes]$preferred_name.2)))
      })
      
    }
  }
  
  if(is.list(modules)){
    
    return(lapply(modules,function(module)PermutString(gene=gene,
                                                       modules=module,
                                                       ppi=ppif,
                                                       connex_genes=connex_genes,
                                                       combined_score.thr=combined_score.thr,
                                                       n.perm=n.perm,
                                                       pct.gap.nconnex=pct.gap.nconnex,
                                                       addNeigbConnex=addNeigbConnex)))
  }else{
    module=modules
  }
  
  #get the number of genes linked to the gene
  n_connex_tot=length(connex_genes[[gene]])
  message('n connexion total for ',gene,': ',n_connex_tot)
  
  #get similar genes
  n_connex_genes<-sapply(connex_genes,length)
  similar_genes<-setdiff(names(which(n_connex_genes>n_connex_tot-n_connex_tot*pct.gap.nconnex&n_connex_genes<n_connex_tot+n_connex_tot*pct.gap.nconnex)),c(gene,module))
  
  message(length(similar_genes),'similar genes connectivity found')
  sampled_genes<-sample(similar_genes,size =min(c(n.perm,length(similar_genes))),replace = F)
  n_connex_sampled_genes<-n_connex_genes[sampled_genes]
  
  #get number of gene of int connected to these genes
  n_connex_obs<-length(intersect(connex_genes[[gene]],module))
  message('n connexion with genes of int: ',n_connex_obs)
  
  n_connex_rand<-sapply(connex_genes[sampled_genes],function(x)length(intersect(x,module)))
  message('n connexion for similar genes: ',paste(head(n_connex_rand,100),collapse = ' '))
  
  return(list(pct.con.permuts=n_connex_rand/n_connex_sampled_genes,
              permuts_genes=sampled_genes,
              pval=sum(n_connex_rand/n_connex_sampled_genes>=n_connex_obs/n_connex_tot)/n.perm,
              gene=gene,
              pct.con.gene=n_connex_obs/n_connex_tot))
  
}
