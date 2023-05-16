

#' Title
#'
#' @param filepath
#' @param bin
#'
#' @return
#' @export
#'
#' @importFrom Biostrings readDNAStringSet subseq letterFrequency
#' @examples
binGC<-function(filepath,bin=10000){
  if(class(bin)!="numeric") { stop("Input bin should be numeric") }
  try<-try(readDNAStringSet(filepath,format="fasta"))
  if ("try-error" %in% class(try)){
    stop("Check the format of your input fasta file !!!")
  }else{
    DNA<-readDNAStringSet(filepath,format="fasta")
  }

  DNAdf<-as.data.frame(DNA@ranges)
  DNAdf<-DNAdf[c("start","end","width","names")]

  info_c<-c()
  for (i in 1:nrow(DNAdf)){

    start=1
    while(start<DNAdf[i,3]){
      end=start+bin-1
      if (end <= DNAdf[i,3]){
        end=end
      }else if(end>DNAdf[i,3]){
        end=DNAdf[i,3]
      }
      info_row<-c(DNAdf[i,4],start,end,i)
      info_c<-rbind(info_c,info_row)
      start=start+bin

    }
  }


  bin_gc_c<-c()
  for (i in 1:nrow(info_c)){
    chr_position<-as.integer(info_c[i,4])
    chr<-info_c[i,1]
    start_p=as.integer(info_c[i,2])
    end_p=as.integer(info_c[i,3])
    gcseq<-subseq(DNA[chr_position],start_p,end_p)
    bin_gc<-letterFrequency(gcseq, "GC", as.prob = TRUE)
    bin_gc_c<-rbind(bin_gc_c,c(chr,start_p,end_p,bin_gc))
  }

  colnames(bin_gc_c)<-c("Name","Start","End","GC")
  bin_gc_c<-as.data.frame(bin_gc_c)
  return(bin_gc_c)
}

