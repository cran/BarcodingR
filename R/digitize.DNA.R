
#' Digitize DNAbin
#' 
#' @description Digitize an object of DNAbin.
#' @param  seqs an object of DNAbin.
#' @return a numeric matrix of DNA sequences digitized.
#' @keywords digitize.DNA
#' @export
#' 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA.
#' @references zhangab2008(at)mail.cnu.edu.cn
#' @examples
#' 
#' data(TibetanMoth)
#' digitized.DNA<-digitize.DNA(seqs=TibetanMoth)
#' digitized.DNA
#' 
#' 


######################################################
########### function: digitalize.DNA start ############
######################################################
digitize.DNA<-function(seqs){
  
  locus<-toupper(as.character(seqs))
  digitized.DNA<-locus
  digitized.DNA[digitized.DNA=="A"]<-1
  digitized.DNA[digitized.DNA=="T"]<-2
  digitized.DNA[digitized.DNA=="G"]<-3
  digitized.DNA[digitized.DNA=="C"]<-4
  digitized.DNA[digitized.DNA=="-"]<-5
  digitized.DNA[digitized.DNA=="N"]<-6
  digitized.DNA[digitized.DNA=="R"]<-0
  digitized.DNA[digitized.DNA=="Y"]<-0
  digitized.DNA[digitized.DNA=="M"]<-0
  digitized.DNA[digitized.DNA=="K"]<-0
  digitized.DNA[digitized.DNA=="S"]<-0
  digitized.DNA[digitized.DNA=="W"]<-0
  digitized.DNA[digitized.DNA=="H"]<-0
  digitized.DNA[digitized.DNA=="B"]<-0
  digitized.DNA[digitized.DNA=="V"]<-0
  digitized.DNA[digitized.DNA=="D"]<-0
  
  
  digitized.DNA<-as.numeric(digitized.DNA)
  #digitized.DNA<-as.matrix(digitized.DNA)
  digitized.DNA2<-array(digitized.DNA,dim=dim(seqs))
  dim(digitized.DNA2)
  
  return(digitized.DNA2)
}


######################################################
########### function: digitalize.DNA end ############
######################################################