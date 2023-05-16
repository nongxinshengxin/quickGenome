


#' Title
#'
#' @param filepath
#' @param format
#' @param telo
#' @param threshold
#'
#' @return
#' @export
#'
#' @importFrom Biostrings readDNAStringSet toString reverse
#' @importFrom stringr str_locate_all str_detect
#' @examples
find_telo <- function(filepath,telo="CCCTAA",threshold=1000,minRepeatNum=2) {
  if(class(telo)!="character") { stop("Input telo should be character") }
  if(class(threshold)!="numeric") { stop("Input threshold should be numeric") }
  if ("try-error" %in% class(try)){
    stop("Check the format of your input fasta file !!!")
  }else{
    DNA<-readDNAStringSet(filepath,format="fasta")
  }
  A2TandG2C<-c("A"="T","G"="C","C"="G","T"="A","N"="N")
  telo_f<-paste(rep(telo,minRepeatNum),collapse  = "")
  telo_r<-paste(A2TandG2C[unlist(strsplit(reverse(telo_f),""))],collapse = "")
  DNAdf<-as.data.frame(DNA@ranges)
  DNAdf<-DNAdf[c("start","end","width","names")]
  left_vector<-c()
  right_vector<-c()
  chr_vector<-c()
  len_vector<-c()
  for (i in 1:length(DNA)){
    all_locate<-c(unlist(str_locate_all(toString(DNA[i]),telo_f)),unlist(str_locate_all(toString(DNA[i]),telo_r)))
    min_p<-min(all_locate)
    max_p<-max(all_locate)
    left<-threshold
    right<-DNAdf[i,3]-threshold
    if (min_p < left){
      min_p=min_p
    }else{
      min_p<-NA
    }
    if (max_p > right){
      max_p=max_p
    }else{
      max_p<-NA
    }
    left_vector<-append(left_vector,min_p)
    right_vector<-append(right_vector,max_p)
    chr_vector<-append(chr_vector,DNAdf[i,4])
    len_vector<-append(len_vector,DNAdf[i,3])
  }
  telo_position<-data.frame(chr=chr_vector,length=len_vector,left=left_vector,right=right_vector)
  return(telo_position)
}
