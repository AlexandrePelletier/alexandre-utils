
#from gencode gtf, get the ensembl id to gene_name links
source('../../utils/r_utils.R')

gtf<-fread('/projectnb/tcwlab/RefData/gencode/hg38/gencode.v45.annotation.gtf')
gtf<-gtf[str_detect(V9,'ENSG')][V3=='gene']
gtf[,gene_id:=str_extract(V9,'ENSG[0-9]+')]
gtf[,gene_name:=str_extract(V9,'gene_name ["A-Za-z0-9-]+')]
gtf[,gene_name:=sapply(gene_name,function(x)strsplit(x,'"')[[1]][2])]
gtf
fwrite(gtf[,.(gene_id,gene_name)],
       '/projectnb/tcwlab/RefData/gencode/hg38/gencode.v45.gene_id_names.tsv')

#mouse
system('wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz')

gtf<-fread('/projectnb/tcwlab/RefData/gencode/mouse/gencode.vM35.annotation.gtf.gz')
gtf<-gtf[str_detect(V9,'ENSM')][V3=='gene']
gtf[,gene_id:=str_extract(V9,'ENSMUSG[0-9]+')]
gtf[,gene_name:=str_extract(V9,'gene_name ["A-Za-z0-9-]+')]
gtf[,gene_name:=sapply(gene_name,function(x)strsplit(x,'"')[[1]][2])]
gtf
fwrite(gtf[,.(gene_id,gene_name)],
       '/projectnb/tcwlab/RefData/gencode/mouse/gencode.vM35.gene_id_names.tsv',sep='\t')
