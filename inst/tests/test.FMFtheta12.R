#' calculate intraspecific maximumal genetic distance theta1 and minimal interspecific distance (theta2) for 
#' a reference dataset with taxon information (species names/IDs)
#' 
#' @param  seqsRef an object of DNAbin in m2.fas format (>seqnames1,species names) or an object of kohonen (som)
#' @return a data.frame with four items, PS(potential species), NN(nearest neighbor), corresponding theta1, theta2, and popSize.PS
#' @keywords FMFtheta12
#' @export 
#' @import adegenet
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
#' http://onlinelibrary.wiley.com/journal/10.1111/(ISSN)1365-294X/earlyview
#' @note  the calculation of theta1 and theta2 is based on euclidan distance ,instead of K2P in this function.
#' @examples
#' require(adegenet)
#' data(seqsSimu) 
#' FMF.theta12<-FMFtheta12(seqsSimu)
#' FMF.theta12
#' 
#' 
library(testthat)
context("FMFtheta12: calculate intraspecific maximumal genetic distance theta1 and minimal interspecific distance")
test_that("FMFtheta12: calculate intraspecific maximumal genetic distance theta1 and minimal interspecific distance",{
  #require(adegenet)
  #setwd("C:/R/myRprojects/ABZ_packages/SOAIC/R")
  #seqs0 <-fasta2DNAbin("locus1.fas") ### important! 
  data(seqsSimu)
  expect_that(class(FMFtheta12(seqsSimu)),
              equals("data.frame")) 
  
})
