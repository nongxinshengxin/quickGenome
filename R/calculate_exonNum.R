

#' Title
#'
#' @param filepath
#' @param format
#' @param collapse
#'
#' @return
#' @export
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @examples
calculate_exonNum<-function(filepath,format="gff",collapse=FALSE){
  if(!format %in% c("gff", "GFF", "gtf", "GTF")) { stop("Format argument needs to be GFF or GTF") }
  gffdb <- makeTxDbFromGFF(filepath, format = format)
  exon.list <- exonsBy(gffdb, by = "tx")
  exon.num <- unlist(lapply(exon.list, FUN = function(x) {NROW(x)}))
  exon.tbl <- table(exon.num)
  if (!collapse) {
    exonNum.df <- data.frame(exon.tbl)
  }else{
    exon.tbl<-as.data.frame(exon.tbl)
    exon.tbl$exon.num<-as.vector(as.integer(exon.tbl$exon.num))
    exonMoreThan10 <- sum(exon.tbl[which(exon.tbl$exon.num>10),2])
    exonNum.df <- rbind(exon.tbl[which(exon.tbl$exon.num<=10),],c(">10",exonMoreThan10))
  }
  return(exonNum.df)

}

