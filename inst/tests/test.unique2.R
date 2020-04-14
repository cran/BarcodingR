#' returns new defined names for unique row
#' @param x a numerical matrix with rownames
#' @param namesBySOM a numberic vector indicating species membership via som for input seqs with unique haplotype,i.e., unique(x or digitized x) has also been used for input of som or others.
#' @return a 2*n matrix, the 1st column the original labels/IDs/rownames of the matrix,  the second unique ID/names/lables
#' @keywords unique2
#' @export 
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created on 2014/7/30 15:50:21
#' @description  unique2() returns a 2-column character matrix where the first corresponds the original labels/IDs of the
#' input matrix (rownames, the input matrix must contain rownames!), the second unique ID/names/lables returned
#' by the function unique() with default parameter settings, or new IDs by som or other methods
#' @details returns a 2-column character matrix where the first corresponds the original labels/IDs,the second unique ID/names/lables returned by the function unique() with default parameter settings, or new IDs by som or other methods
#' 
#' @examples
#' x<-array(1:12,c(3,4))
#' rownames(x)<-c("s1","s2","s3")
#' x<-rbind(x,x[3,])
#' rownames(x)[4]<-"s4"
#' namesBySOM<-c(1,2,2)
#' out<-unique2(x)
#' out
#' out<-unique2(x,namesBySOM)
#' out


library(testthat)
context("unique2: returns a 2-column character matrix where the first corresponds the original labels/IDs,the second unique ID...")
test_that("unique2: returns a 2-column character matrix where the first corresponds",{
  
  x<-array(1:12,c(3,4))
  rownames(x)<-c("s1","s2","s3")
  x<-rbind(x,x[3,])
  rownames(x)[4]<-"s4"
  namesBySOM<-c(1,2,2)
  out<-unique2(x,namesBySOM)
  expect_that(class(out),
              equals("matrix")) 
})
