
PrepDESeq2Design<-function(mtd,covs,scale_num=FALSE){
  
  #remove if all sample have the same value
  to_keep<-sapply(covs,function(cov)length(unique(mtd[[cov]]))!=1)
  if(sum(!to_keep)>1)message('remove ',paste(covs[!to_keep],collapse = ', '),'because all samples have the same value')
  covs_filtered<-covs[to_keep]
  
  #if a factor, need to have at least 2 sample in each category
  to_keep<-sapply(covs_filtered,function(cov)ifelse(!is.numeric(mtd[[cov]]),all(table(mtd[[cov]])!=1),TRUE))
  if(sum(!to_keep)>1)message('remove ',paste(covs_filtered[!to_keep],collapse = ', '),' factor(s) because need to have at least 2 samples in each category')
  
  covs_filtered<-covs_filtered[to_keep]
  
  #if numerical, scales the data
  
  if(scale_num){
    numerical_factors<-covs_filtered[mtd[,sapply(.SD,is.numeric),.SDcols=covs_filtered]]
    
    if(length(numerical_factors)>0){
      numerical_factors_scaled<-paste0(numerical_factors,'_scaled')
      message('scaling ',paste(numerical_factors,collapse=', '),' in new column ',paste(numerical_factors_scaled,collapse=', '))
      mtd[,(numerical_factors_scaled):=lapply(.SD,scale),.SDcols=numerical_factors]
      covs_filtered<-sapply(covs_filtered, function(x)ifelse(x%in%numerical_factors,paste0(x,'_scaled'),x))
    }
    
  }
 
  formul_clean<-formula(paste0("~", paste(covs_filtered,collapse ='+')))
  return(formul_clean)

  
}

ScaleDESeq2Covs<-function(mtd,covs){
  mtd_scaled<-copy(mtd)
  numerical_factors<-covs[mtd_scaled[,sapply(.SD,is.numeric),.SDcols=covs]]
  numerical_factors_scaled<-paste0(numerical_factors,'_scaled')
  message('scaling ',paste(numerical_factors,collapse=', '),' in new column ',paste(numerical_factors_scaled,collapse=', '))
  mtd_scaled[,(numerical_factors_scaled):=lapply(.SD,scale),.SDcols=numerical_factors]
  return(mtd_scaled)
}
