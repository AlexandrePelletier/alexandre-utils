#seurat utils
require('Seurat')

splitCluster<-function(sample,cluster,dimensions=1:6,resol=0.3){
  library(stringr)
  cluster_to_div<-subset(sample,idents = cluster)
  
  cluster_to_div<-RunPCA(cluster_to_div,features=rownames(cluster_to_div))
  cluster_to_div<-RunUMAP(cluster_to_div,dims=dimensions)
  cluster_to_div<-FindNeighbors(cluster_to_div,dims = dimensions)
  cluster_to_div<-FindClusters(cluster_to_div,algorithm = 1,resolution = resol)
  Idents(sample,cells=colnames(cluster_to_div))<-paste(cluster,Idents(cluster_to_div),sep=".")
  
  
  return(sample)
}


DownSampleSeurat<-function(object,prop,group.by=NULL,min.cells=20){
  require('data.table')
  
  mtd<-data.table(object@meta.data,keep.rownames = 'bc')
  mtd[,to_keep:=bc%in%sample(bc,size = ifelse(prop*.N>min.cells,round(prop*.N),min.cells)),by=group.by]
  object_sub<-object[,mtd[(to_keep)]$bc]
  
  return(object_sub)
}
