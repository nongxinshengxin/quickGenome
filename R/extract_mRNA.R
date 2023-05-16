

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
extract_mRNA<-function(fasta,gff,outputfile){
  options(warn = -1)
  DNA<-readDNAStringSet(fasta,format="fasta")

  gffdf<-read_delim(gff,delim = "\t",col_names = F)%>%
    drop_na()%>%
    filter_all(any_vars(!str_detect(., pattern = "#")))%>%
    select(c(1,3,4,5,7,9))%>%
    rename_with(~c("chr","type","start","end","strand","attribution"),c(1:6))


  if (!"gene" %in% gffdf$type | !"CDS" %in% gffdf$type | !"exon" %in% gffdf$type){
    stop("GFF format is wrong, please check the format!")
  }


  exondf<-filter(gffdf,type=="exon")
  exondf<-exondf%>%mutate(idname=gsub("Parent=","",grep("Parent.*",unlist(strsplit(exondf$attribution,";")),value = T)))

  exon_c<-c()
  for (i in 1:nrow(exondf)){
    chr=as.character(exondf[i,1])
    start=as.integer(exondf[i,3])
    end=as.integer(exondf[i,4])
    strand=exondf[i,5]
    id=as.character(exondf[i,7])
    seq=subseq(DNA[chr],start=start,end=end)
    seq=toString(seq)
    exon_c<-rbind(exon_c,c(id,seq,strand))
  }

  exon_c<-as_tibble(exon_c)%>%pivot_wider(names_from = V1,values_from = V2,values_fill = NA)

  id_c<-c()
  seq_c<-c()
  if (exon_c[1,1]=="+"){
    for (i in 2:ncol(exon_c)){
      id<-colnames(exon_c[,i])
      id_c<-append(id_c,id)
      if (!is.null(unlist(exon_c[1,i]))){
        seq=paste(unlist(exon_c[1,i]),collapse = "")
      }
      if (!is.null(unlist(exon_c[2,i]))){
        seq=paste(unlist(exon_c[2,i]),collapse = "")
        seq=DNAString(seq)
        seq=reverseComplement(seq)
        seq=toString(seq)
      }
      seq_c<-append(seq_c,seq)

    }
  }else{
    for (i in 2:ncol(exon_c)){
      id<-colnames(exon_c[,i])
      id_c<-append(id_c,id)
      if (!is.null(unlist(exon_c[2,i]))){
        seq=paste(unlist(exon_c[2,i]),collapse = "")
      }
      if (!is.null(unlist(exon_c[1,i]))){
        seq=paste(unlist(exon_c[1,i]),collapse = "")
        seq=DNAString(seq)
        seq=reverseComplement(seq)
        seq=toString(seq)
      }
      seq_c<-append(seq_c,seq)
      # seq_c<-append(seq_c,seq2)
    }
  }

  res<-paste0(">",id_c,"\n",seq_c)
  write.table(res,file = outputfile,row.names = F,quote = F, col.names=F)


}

