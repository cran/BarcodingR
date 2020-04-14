
#' Calculate euclidean distance of two numeric vectors
#' 
#' @param  v1 a numeric vector
#' @param  v2 another numeric vector
#' @return one numeric
#' @keywords eucl.dist.two.vect
#' @export 
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created on 2014/4/20 15:11:03 by ABZ,CNU,CHINA
#' @examples
#' v1<-1:5
#' v2<-6:10
#' eucl.dist.two.vect(v1,v2)
#' 
#' 


library(testthat)
context("eucl.dist.two.vect: calculate euclidean distance of two numeric vectors")
test_that("eucl.dist.two.vect: calculate euclidean distance of two numeric vectors",{
  v1<-1:5
  v2<-6:10
  #expect_that(class(DNAbin2kmerFreqMatrix(seqs=woodmouse,kmer=3)),equals("matrix")) 
  expect_that(eucl.dist.two.vect(v1,v2),equals(11.18034)) 
})

