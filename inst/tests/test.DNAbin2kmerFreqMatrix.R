

#' Convert an object of DNAbin into kmer frequence matrix by specifying a certain length of kmer'
#' 
#' @param  seqs an object of DNAbin (may not be aligned!)
#' @param  kmer an integer
#' @return a matrix of kmer frequency
#' @keywords DNAbin2kmerFreqMatrix
#' @export 
#' @import seqinr
#' @import ape
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created on on 2014/8/20 17:19:23
#' @examples
#' 
#' require(seqinr)
#' require(ape)
#' data(woodmouse)
#' kmer.Freq<-DNAbin2kmerFreqMatrix(seqs=woodmouse,kmer=3)
#' kmer.Freq
#' 
#' 



library(testthat)
context("DNAbin2kmerFreqMatrix: convert object of DNAbin into kmer frequence matrix")
test_that("DNAbin2kmerFreqMatrix: convert object of DNAbin into kmer frequence matrix",{
  require(seqinr)
  require(ape)
  data(woodmouse)
  expect_that(class(DNAbin2kmerFreqMatrix(seqs=woodmouse,kmer=3)),equals("matrix")) 
  
})



