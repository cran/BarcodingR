#' calculate fuzzy membership function (FMF) value
#' 
#' @param  xtheta12 a numeric vector of 1*3, with the 1st element (x: the genetic distance of an query
#' DNA sequence to a species, the 2nd, theta1 - the maximum intraspecific genetic distance of the 
#' the potential species assigned (PS), the 3rd, theta2 - the interspecific species between the  
#' PS and its nearest neighbor (NN).
#' @return fuzzy membership function (FMF) value in the range of [0,1], 0 indicating not belong to the PS
#' 1 indicating completely belonging PS, other values indicating the degree of belonging to PS.
#' @keywords FMF
#' @export 
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
#' http://onlinelibrary.wiley.com/journal/10.1111/(ISSN)1365-294X/earlyview
#' @note  initially created on 2014/8/23 10:57:59
#' @examples
#' xtheta12<-c(0.8,0.7,1.2)
#' FMF<-FMF(xtheta12)
#' FMF


library(testthat)
context("FMF: calculate fuzzy membership function (FMF) value")
test_that("FMF: calculate fuzzy membership function (FMF) value",{
  xtheta12<-c(0.8,0.7,1.2)
  expect_that(FMF(xtheta12),
              equals(0.92)) 
  
})

