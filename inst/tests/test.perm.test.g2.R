
#' Check whether significant differentiation occurred between two groups (matrices)
#' 
#' @param  g1 a numeric matrix with row indicating different specimens, columns indicating different variables
#' @param  g2 a numeric matrix with row indicating different specimens, columns indicating different variables
#' @param  dist.method a character string, with default value of "euclidian"
#' @param  group.method a character string, with default value of "centroid"
#' @param  nRep the number of permutation, with default value of 99
#' @return a two-elements vector indicating significance (1-significant,0-not significant) and P-value
#' @keywords perm.test.g2
#' @export 
#' @import ape
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created 2014/8/12 14:59:50
#' @description Check whether significant differentiation occurred between two groups (matrices) with permutation test.
#' @examples
#' 
#' require(ape)
#' data(woodmouse)
#' digitized.DNA<-digitize.DNA(seqs=woodmouse)
#' g1<-digitized.DNA[1:5,]
#' g2<-digitized.DNA[6:10,]
#' out.perm.t<-perm.test.g2(g1,g2,dist.method = "euclidian", group.method = "centroid",nRep = 99)
#' out.perm.t




library(testthat)
context("perm.test.g2: Check whether significant differentiation occurred between two groups (matrices)")
test_that("perm.test.g2: Check whether significant differentiation occurred between two groups (matrices)",{
  require(ape)
  data(woodmouse)
  digitized.DNA<-digitize.DNA(seqs=woodmouse)
  g1<-digitized.DNA[1:5,]
  g2<-digitized.DNA[6:10,]
  expect_that(length(perm.test.g2(g1,g2,dist.method = "euclidian", 
                           group.method = "centroid",nRep = 99)),
              equals(2)) 
  
})