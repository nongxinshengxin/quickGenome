
#' Title
#'
#' @param fasta
#' @param genelist
#' @param outputfile
#' @param Negate
#' @param Str
#'
#' @return
#' @export
#'
#' @examples
extract_fasta<-function(fasta,genelist,outputfile,Negate=FALSE,Str=FALSE){
  DNA<-readDNAStringSet(fasta,format="fasta")
  DNAdf<-as.data.frame(DNA@ranges)
  DNAdf<-DNAdf%>%select(start,end,width,names)%>%
    mutate(position=str_split(DNAdf$names, " ", simplify = TRUE)[,1])

  if (!Str){
    listdf<-read_delim(genelist,delim = "\t",col_names = F)
  }else{
    listdf<-as_tibble(genelist)
  }


  id_c<-c()
  seq_c<-c()
  if (!Negate){
    for (i in 1:nrow(listdf)) {
      idname=as.character(listdf[i,1])
      chr=DNAdf[which(DNAdf$position==idname),4]
      if (length(chr) ==0){
        warning(paste(idname,"is not matching",sep = " "))
      }else{
        id_c<-append(id_c,idname)
        seq=toString(DNA[chr])
        seq_c<-append(seq_c,seq)
      }
    }
  }else{
    for (i in 1:nrow(listdf)) {
      idname=as.character(listdf[i,1])
      chr=DNAdf[which(DNAdf$position==idname),4]
      if (length(chr) ==0){
        warning(paste(idname,"is not matching",sep = " "))
      }else{
        id_c<-append(id_c,chr)
      }
    }
    id_c<-setdiff(DNAdf$names,id_c)
    for (i in 1:length(id_c)){
      seq=toString(DNA[id_c[i]])
      seq_c<-append(seq_c,seq)
    }
  }


  res<-paste0(">",id_c,"\n",seq_c)
  write.table(res,file = outputfile,row.names = F,quote = F, col.names=F)
}

