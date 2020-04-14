

#' Divide an integer into two integers,which product is almost equal to the input integer,and with minal difference between these two numbers
#' 
#' @param  one an integer
#' @return a vector with two elements of integer
#' @keywords one2Two
#' @export 
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created on on 2014/8/20 17:19:23
#' @examples
#' one2Two(45)
#' 


library(testthat)
context("one2Two: Divide an integer into two integers")
test_that("one2Two: Divide an integer into two integers",{
  expect_that(one2Two(45),equals(c(5,9))) 
})





