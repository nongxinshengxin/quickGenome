#' Title
#'
#' @param fasta
#' @param gff
#' @param outputfile
#'
#' @return
#' @export
#'
#' @examples
extract_gene<-function(fasta,gff,outputfile){
  DNA<-readDNAStringSet(fasta,format="fasta")
  DNAdf<-as.data.frame(DNA@ranges)
  DNAdf<-DNAdf%>%select(start,end,width,names)%>%
    mutate(position=str_split(DNAdf$names, " ", simplify = TRUE)[,1])

  gffdf<-read_delim(gff,delim = "\t",col_names = F)%>%
    drop_na()%>%
    filter_all(any_vars(!str_detect(., pattern = "#")))%>%
    select(c(1,3,4,5,7,9))%>%
    rename_with(~c("chr","type","start","end","strand","attribution"),c(1:6))

  if (!"gene" %in% gffdf$type){
    stop("GFF format is wrong, please check the format!")
  }

  genedf<-filter(gffdf,type=="gene")
  genedf<-genedf%>%mutate(idname=gsub("ID=","",grep("ID.*",unlist(strsplit(genedf$attribution,";")),value = T)))


  id_c<-c()
  seq_c<-c()
  for (i in 1:nrow(genedf)){
    chr=as.character(genedf[i,1])
    chr=DNAdf[which(DNAdf$position==chr),4]
    start=as.integer(genedf[i,3])
    end=as.integer(genedf[i,4])
    strand=genedf[i,5]
    id=as.character(genedf[i,7])
    id_c<-append(id_c,id)
    seq=subseq(DNA[chr],start=start,end=end)
    if (strand=="+"){
      seq=toString(seq)
    }else{
      seq=reverseComplement(seq)
      seq=toString(seq)
    }
    seq_c<-append(seq_c,seq)
  }

  res<-paste0(">",id_c,"\n",seq_c)
  write.table(res,file = outputfile,row.names = F,quote = F, col.names=F)

}

