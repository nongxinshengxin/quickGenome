
library(tidyverse)
library(Biostrings)


#' Title
#'
#' @param fasta
#' @param gff
#' @param outputfile
#' @param translation
#'
#' @return
#' @export
#'
#' @importFrom readr read_delim
#' @importFrom Biostrings readDNAStringSet toString reverse subseq DNAString translate
#' @importFrom tibble tibble
#' @importFrom tidyr drop_na pivot_wider
#' @importFrom dplyr `%>%` filter select rename_with mutate any_vars filter_all
#' @examples
extract_CDS<-function(fasta,gff,outputfile,translation=FALSE,if.fuzzy.codon="error"){
  options(warn = -1)
  DNA<-readDNAStringSet(fasta,format="fasta")
  DNAdf<-as.data.frame(DNA@ranges)
  DNAdf<-DNAdf%>%select(start,end,width,names)%>%
    mutate(position=str_split(DNAdf$names, " ", simplify = TRUE)[,1])

  gffdf<-read_delim(gff,delim = "\t",col_names = F)%>%
    drop_na()%>%
    filter_all(any_vars(!str_detect(., pattern = "#")))%>%
    select(c(1,3,4,5,7,9))%>%
    rename_with(~c("chr","type","start","end","strand","attribution"),c(1:6))


  if (!"CDS" %in% gffdf$type){
    stop("GFF format is wrong, please check the format!")
  }


  cdsdf<-filter(gffdf,type=="CDS")
  cdsdf<-cdsdf%>%mutate(idname=gsub("Parent=","",grep("Parent.*",unlist(strsplit(cdsdf$attribution,";")),value = T)))

###### R的for循环太慢了
  # cds_c<-c()
  # for (i in 1:nrow(cdsdf)){
  #   chr=as.character(cdsdf[i,1])
  #   chr=DNAdf[which(DNAdf$position==chr),4]
  #   start=as.integer(cdsdf[i,3])
  #   end=as.integer(cdsdf[i,4])
  #   strand=cdsdf[i,5]
  #   id=as.character(cdsdf[i,7])
  #   seq=subseq(DNA[chr],start=start,end=end)
  #   seq=toString(seq)
  #   cds_c<-rbind(cds_c,c(id,seq,strand))
  # }
  #
  # cds_c<-as_tibble(cds_c)%>%pivot_wider(names_from = V1,values_from = V2,values_fill = NA)

  func <- function(x){
    chr=as.character(x[1])
    chr=DNAdf[which(DNAdf$position==chr),4]
    start=as.integer(x[3])
    end=as.integer(x[4])
    strand=x[5]
    id=as.character(x[7])
    seq=subseq(DNA[chr],start=start,end=end)
    seq=toString(seq)
    cds_c<-c(id,seq,strand)
    return(cds_c)
  }


  a3 <- apply(cdsdf,1,func)
  a3<-t(a3)
  cds_c<-as_tibble(a3)%>%pivot_wider(names_from = V1,values_from = V2,values_fill = NA)

  id_c<-c()
  seq_c<-c()
  if (!translation){
    if (cds_c[1,1]=="+"){
      for (i in 2:ncol(cds_c)){
        id<-colnames(cds_c[,i])
        id_c<-append(id_c,id)
        if (!is.null(unlist(cds_c[1,i]))){
          seq=paste(unlist(cds_c[1,i]),collapse = "")
        }
        if (!is.null(unlist(cds_c[2,i]))){
          seq=paste(unlist(cds_c[2,i]),collapse = "")
          seq=DNAString(seq)
          seq=reverseComplement(seq)
          seq=toString(seq)
        }
        seq_c<-append(seq_c,seq)

      }
    }else{
      for (i in 2:ncol(cds_c)){
        id<-colnames(cds_c[,i])
        id_c<-append(id_c,id)
        if (!is.null(unlist(cds_c[2,i]))){
          seq=paste(unlist(cds_c[2,i]),collapse = "")
        }
        if (!is.null(unlist(cds_c[1,i]))){
          seq=paste(unlist(cds_c[1,i]),collapse = "")
          seq=DNAString(seq)
          seq=reverseComplement(seq)
          seq=toString(seq)
        }
        seq_c<-append(seq_c,seq)
        # seq_c<-append(seq_c,seq2)
      }
    }
  }else{
    if (cds_c[1,1]=="+"){
      for (i in 2:ncol(cds_c)){
        id<-colnames(cds_c[,i])
        id_c<-append(id_c,id)
        if (!is.null(unlist(cds_c[1,i]))){
          seq=paste(unlist(cds_c[1,i]),collapse = "")
          seq=DNAString(seq)
          seq=translate(seq,if.fuzzy.codon=if.fuzzy.codon)
          seq=toString(seq)
        }
        if (!is.null(unlist(cds_c[2,i]))){
          seq=paste(unlist(cds_c[2,i]),collapse = "")
          seq=DNAString(seq)
          seq=reverseComplement(seq)
          seq=translate(seq,if.fuzzy.codon=if.fuzzy.codon)
          seq=toString(seq)
        }
        seq_c<-append(seq_c,seq)

      }
    }else{
      for (i in 2:ncol(cds_c)){
        id<-colnames(cds_c[,i])
        id_c<-append(id_c,id)
        if (!is.null(unlist(cds_c[2,i]))){
          seq=paste(unlist(cds_c[2,i]),collapse = "")
          seq=DNAString(seq)
          seq=translate(seq,if.fuzzy.codon=if.fuzzy.codon)
          seq=toString(seq)
        }
        if (!is.null(unlist(cds_c[1,i]))){
          seq=paste(unlist(cds_c[1,i]),collapse = "")
          seq=DNAString(seq)
          seq=reverseComplement(seq)
          seq=translate(seq,if.fuzzy.codon=if.fuzzy.codon)
          seq=toString(seq)
        }
        seq_c<-append(seq_c,seq)
        # seq_c<-append(seq_c,seq2)
      }
    }
  }




  res<-paste0(">",id_c,"\n",seq_c)
  write.table(res,file = outputfile,row.names = F,quote = F, col.names=F)

}


