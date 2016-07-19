
#' Extracts Lables of Samples
#' 
#' @description Extract sequence names from different objects of DNAbin, including generated from fasta2DNAbin() (package:adegenet),
#' and read.dna() (package:ape).
#' 
#' @param  seqs object of class "DNAbin", generated from fasta2DNAbin() (package:adegenet), and read.dna() (package:ape). 
#' @return a character string array/vector.
#' @keywords NAMES
#' @export 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA. 
#' @references zhangab2008(at)mail.cnu.edu.cn
#' @examples
#' 
#' data(TibetanMoth) 
#' seqNames<-NAMES(TibetanMoth)
#' seqNames


NAMES<-function(seqs){
  
  if(mode(seqs)=="raw"){
  SeqNames<-attr(seqs,"dimnames")[[1]]
  }else{ ### mode(seqs)=="list"
    SeqNames<-names(seqs)  
  }
  
  #names(SeqNames)<-NULL
  
  if(length(SeqNames)==0) {
    stop("the mode(seqs) is wrong!")
  }else{
    return(SeqNames) }
}




















