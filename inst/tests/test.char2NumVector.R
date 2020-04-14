
#' Convert an array of characters into a vector of numeric in increasing order.
#' 
#' @param c a character string
#' @return a vector of numeric
#' @keywords char2NumVector
#' @export 
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at
#' zhangab2008@mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created on 2014/8/1 8:42:55
#' @examples
#' char2NumVector()
#' c<-c("A","B","Z")
#' n<-char2NumVector(c)
#' 
#' 
#' 

library(testthat)
context("char2NumVector: input an array of characters")
test_that("char2NumVector: input an array of characters",{
  
  expect_that(char2NumVector(c("A","B","Z")),equals(c(1,2,3))) 
  
})
