

#' Title
#'
#' @param filepath
#' @param format
#'
#' @return
#' @export
#'
#' @importFrom Biostrings readDNAStringSet letterFrequency
#' @examples
genome_basic_Info<-function(filepath){
  try<-try(readDNAStringSet(filepath,format="fasta"))
  if ("try-error" %in% class(try)){
    stop("Check the format of your input fasta file !!!")
  }else{
    DNA<-readDNAStringSet(filepath,format="fasta")
  }
  DNAdf<-as.data.frame(DNA@ranges)
  DNAdf<-DNAdf[c("start","end","width","names")]
  DNAdf<-DNAdf[order(DNAdf$width,decreasing = T),]


  sumlen<-sum(DNAdf$width)
  print(paste("The size of genome is",sumlen))
  contignum<-nrow(DNAdf)
  print(paste("The number of contig/scaffold/chromosome is",contignum))
  Largest_contig<-DNAdf[1,3]
  print(paste("The size of largest contig/scaffold/chromosome is",Largest_contig))
  GC<-sum(letterFrequency(DNA, "GC"))/sumlen
  print(paste("GC content is",GC))
  lenplus=0
  for (i in 1:nrow(DNAdf)){
    lenplus=lenplus+DNAdf[i,3]
    ratio=lenplus/sumlen
    if (ratio>=0.5){
      print(paste("N50 is",DNAdf[i,3]))
      print(paste("L50 is",i))
      break
    }
  }
  lenplus2=0
  for (i in 1:nrow(DNAdf)){
    lenplus2=lenplus2+DNAdf[i,3]
    ratio=lenplus2/sumlen
    if (ratio>=0.9){
      print(paste("N90 is",DNAdf[i,3]))
      print(paste("L90 is",i))
      break
    }
  }
}
