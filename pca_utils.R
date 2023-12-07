#ALL PCA BASED ANALYSIS UTILITIES


RunPca<-function(norm_mat,scale=TRUE,remove_incomplete_rows=T){
  #return  the PCA of the transposed matrix : if the samples are the columns, factors will be estimated by samples
  if(remove_incomplete_rows){
    norm_mat<-norm_mat[rowSums(is.na(norm_mat))==0,] #remove feature with NA
    
  }
    norm_mat<-norm_mat[which(apply(norm_matf, 1, var)!=0),] #remove feature without variance
    message('removing ',nrow(norm_mat)-nrow(norm_matf),' features (',nrow(norm_matf), 'remaining) because missing value or without variance')
    return(prcomp(t(norm_matf),scale.=scale))
    
}

pctPC<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-paste0('PC',rngPCs)
  return( pct.varPCs)
}


PcaPlot<-function(pca,mtd,group.by, pc_x='PC1', pc_y='PC2',sample_col='sample_id',scale=TRUE,return_pcs_mtd=FALSE){
  
  pca_dt<-merge(data.table(pca$x,keep.rownames = sample_col),mtd)
  
  pctpcs<-pctPC(pca)
  
  ps<-lapply(group.by, function(c){
    p<-ggplot(pca_dt)+geom_point(aes_string(x=pc_x,y=pc_y,col=c))+
      labs(x=paste0(pc_x,' (',round(pctpcs[pc_x]*100),'%)'),
           y=paste0(pc_y,' (',round(pctpcs[pc_y]*100),'%)'))+
      theme_bw()
    return(p+ggtitle(c))
  })
  
  
  if(return_pcs_mtd){
    print(wrap_plots(ps))
    return(pca_dt)
  }
  else{
    return(wrap_plots(ps))
  }
}



CorrelCovarPCs<-function(pca,mtd,
                         sample_col='sample',vars_num=NULL,
                         vars_fac=NULL,rngPCs=1:10,res="pval",seuilP=0.1,return=TRUE,plot.results=FALSE){
  require(data.table)
  
  pcs<-data.frame(pca$x)
  names(rngPCs)<-paste0('PC',rngPCs)
  
  #if(is.null(vars_num))vars_num=colnames(mtd)[sapply(mtd, is.numeric)]
  #if(is.null(vars_fac))vars_fac=colnames(mtd)[!sapply(mtd, is.numeric)&colnames(mtd)!=sample_col]
  
  
  if(is.data.table(mtd)){
    mtd[,(vars_fac):=lapply(.SD,as.factor),.SDcols=vars_fac]
    
    mtd<-data.frame(mtd,row.names = sample_col)
  }
  
  if(length(vars_num)>0){
    mtd_num=mtd[rownames(pcs),vars_num]
    mtd_num=t(mtd_num)
    split_num=split(mtd_num,rownames(mtd_num))
    res_num=lapply(split_num,function(x,ret=res){
      FAC1.p<-rep(0,length(rngPCs))
      FAC1.r2<-rep(0,length(rngPCs))
      for (i in rngPCs){
        FAC1<-as.numeric(x)
        FAC1<-lm(pcs[,i]~FAC1)
        FAC1.p[i]<-anova(FAC1)$Pr[1]
        FAC1.r2[i]<-summary(FAC1)$adj.r.squared
      }
      if(ret=="pval"){
        return(FAC1.p)
      }else if(ret=="r2"){
        return(FAC1.r2)
      }
    })
    
    res.num<-do.call(rbind,res_num)
    
  }else{
    res.num<-data.frame()
  }
  
  
  if(length(vars_fac)>0){
    mtd_fac=mtd[rownames(pcs),vars_fac]
    
    mtd_fac<-mtd_fac[,sapply(mtd_fac, function(x)length(levels(x))>1)]
    mtd_fac=t(mtd_fac)
    split_fac=split(mtd_fac,rownames(mtd_fac))
    res_fac=lapply(split_fac,function(x,ret=res){
      FAC1.p<-rep(0,length(rngPCs))
      FAC1.r2<-rep(0,length(rngPCs))
      for (i in rngPCs){
        FAC1<-as.factor(x)
        FAC1<-lm(pcs[,i]~FAC1)
        FAC1.p[i]<- anova(FAC1)$Pr[1]
        FAC1.r2[i]<-summary(FAC1)$adj.r.squared
      }
      if(ret=="pval"){
        return(FAC1.p)
      }else if(ret=="r2"){
        return(FAC1.r2)
      }})
    res.fac<-do.call(rbind,res_fac)
    
  }else{
    res.fac<-data.frame()
    
  }
  
  final_res<-rbind(res.num,res.fac)
  final_res<-data.matrix(final_res)
  
  if(plot.results){
    require(pheatmap)
    resToPlot<-final_res
    if(res=="pval"){
      
      resToPlot[which(resToPlot>seuilP)]<-1 ####here I basicaly put them to 1 if less than 0.1
      resToPlot<--log10(resToPlot)
      breakRes<-c(40,20,10:1, 0.5,0.1)
    }else{
      resToPlot[resToPlot<0]<-0
      resToPlot<-resToPlot
      breakRes<-NA
      
    }
    
    pct.varPCs<-pctPC(pca,rngPCs)*100
    vars<-rownames(resToPlot)
    
    if(plot.results){
      plotPvalsHeatMap(resToPlot[vars,rngPCs],
                       labels_col= paste0(names(rngPCs),"(",round(pct.varPCs[names(rngPCs)],0),"%)"),
                       col_breaks = breakRes)
      
    }
    
 
  }
  
  
  
  if(return){
    colnames(final_res)<-paste0("PC",rngPCs)
    return(final_res) 
  }
  
  
  
}


plotPvalsHeatMap<-function(pvals_mat,main='-log10(Pvalue)',p.thr=0.1,col_breaks=c(20,10:1, 0.5,0.1),
                           legend_breaks=NA,cluster_rows = F,cluster_cols = F,
                           labels_col=NULL,
                           labels_row=NULL,
                           fontsize=10,
                           fontsize_number = 0.8*fontsize){
  require(pheatmap)
  
  pvals_mat[which(pvals_mat>p.thr)]<-1 #### put them to 1 if less than 0.1
  pvals_mat<--log10(pvals_mat)
  
  vars<-rownames(pvals_mat)

  pheatmap(pvals_mat,main = main,
           cluster_rows = cluster_rows,cluster_cols = cluster_cols,
           display_numbers = T,
           fontsize = fontsize,
           fontsize_number = fontsize_number,
           labels_col=labels_col,
           labels_row=labels_row,
           color = colorRampPalette(c("white", "red"))(length(col_breaks)-1),
           breaks = sort(col_breaks),legend_breaks = legend_breaks)
  
}


correl<-function(x,y=NULL,ret="pval",verbose=F){
  if(is.null(y)){
    y=unlist(x[,2],use.names = F)
    x=unlist(x[,1],use.names = F)
  }
  
  if(is.numeric(x)){
    
    if(verbose){
      message("linear modeling ")
    }
    
    res<-lm(x~y)
    if(verbose){
      print(summary(res))
    }
    if(ret=="r2"){
      
      return(summary(res)$adj.r.squared)
    }else if(ret=="pval"){
      return(anova(res)$Pr[1])
      
    }else{
      return(
        list(p=anova(res)$Pr[1],r2=summary(res)$adj.r.squared)
      )
    }
    
  }else if(all(sapply(list(x,y),is.factor))){
    if(verbose){
      message("Chi-squared independance test")
    }
    
    x<-factor(x,levels = unique(x))
    y<-factor(y,levels = unique(y))
    tableF<-table(x,y)
    test<-chisq.test(tableF)
    if(verbose){
      print(test)
    }
    
    return(test$p.value)
  }
  
}
