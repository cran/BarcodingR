
#' digitize an object of DNAbin.
#' 
#' @param  seqs an object of DNAbin (may not be aligned!)
#' @return a numeric matrix of DNA sequences digita
#' @keywords digitize.DNA
#' @export 
#' @import ape
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.In
#' @examples
#' 
#' require(ape)
#' data(woodmouse)
#' digitized.DNA<-digitize.DNA(seqs=woodmouse)
#' digitized.DNA
#' 
#' 

library(testthat)
context("digitize.DNA: digitize an object of DNAbin.")
test_that("digitize.DNA: digitize an object of DNAbin.",{
  require(ape)
  data(woodmouse)
  expect_that(class(digitize.DNA(seqs=woodmouse)),equals("matrix")) 
              
})







