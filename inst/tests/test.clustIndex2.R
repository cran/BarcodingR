#' calculate clustIndex
#' 
#' @param  y Object of class "cclust" returned by a clustering algorithm such as kmeans
#' @param  x Data matrix where columns correspond to variables and rows to observations
#' @param  index a character string 
#' @return a numeric value of cluster index
#' @keywords clustIndex2
#' @export 
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references modified from clustIndex of package cclust
#' @examples
#' 
#' data(outSom) 
#' y<-kohonenClass2cclust(outSom)
#' x<-outSom$data
#' clIndex2<-clustIndex2(y,x,index = "calinski")
#' clIndex2

#' 
library(testthat)
context("clustIndex2: calculate clustIndex")
test_that("clustIndex2: calculate clustIndex",{
  data(outSom) 
  y<-kohonenClass2cclust(outSom)
  x<-outSom$data
  clIndex2<-clustIndex2(y,x,index = "calinski")
  clIndex2<-floor(clIndex2)
  expect_that(class(clIndex2),
              equals("numeric")) 
})


