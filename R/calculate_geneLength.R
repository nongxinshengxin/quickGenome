
#' Title
#'
#' @param filepath 
#' @param format 
#' @param type 
#' @param bin 
#' @param threshold 
#'
#' @return
#' @export
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF genes transcripts
#' @importFrom Biostrings width
#' @examples
calculate_geneLength<-function(filepath,format="gff",type="genes",bin=500,threshold=4000){
  if(!format %in% c("gff", "GFF", "gtf", "GTF")) { stop("Format argument needs to be GFF or GTF") }
  if(!type %in% c("genes", "transcripts")) { stop("Type argument needs to be genes or transcripts") }
  gffdb <- makeTxDbFromGFF(filepath, format = format)
  if (type=="genes"){
    ge.length <- width(genes(gffdb))
  }else if(type=="transcripts"){
    ge.length <- width(transcripts(gffdb))
  }

  aa = table(cut(ge.length, breaks=seq(0,threshold, bin))) 
  
  info_c=c()
  start=1
  while(start<threshold){
    end=start+bin-1
    info_row<-paste(start,end,sep = "-")
    info_c<-append(info_c,info_row)
    start=start+bin
    
  }
  info_c<-append(info_c,paste(">",as.character(threshold)))
  morethan4000 <- length(ge.length[ge.length > threshold])
  geLen.df <- data.frame(range=info_c,number = c(as.vector(aa), morethan4000))
  return(geLen.df)
}