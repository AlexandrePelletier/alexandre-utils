#melt the matrix expression and metatadata to save a single csv file by gene containing samples metadata
out<-'outputs/02-csv_files_tmm_genes'
dir.create(out)

library(data.table)

expr<-fread('data/tmm_data.csv')

mtd<-fread('data/tmm_meta.csv')

#melt the matrix (wide to long transformation)
expr_dt<-melt(expr,id.vars = 'gname',variable.name = 'id',value.name = 'value')

#merge with metadata
expr_dt<-merge(expr_dt,
               mtd[,.(id,SEX,cell_type,APOE,indiv_id,proportion)]) #keep only useful metadata

#split bu genes and save
expr_dt_list<-split(expr_dt,by='gname')

lapply(names(expr_dt_list), function(gene_name){
  fwrite(expr_dt_list[[gene_name]],file.path(out,paste0('tmm_',gene_name,'.csv')))
})
