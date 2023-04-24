#seurat utils
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
