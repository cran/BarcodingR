#' transfers  "kohonen" class into "cclust" class
#' 
#' @param  x Object of class "kohonen" returned by a clustering algorithm such as som
#' @return object of class "cclust"
#' @keywords kohonenClass2cclust
#' @export 
#' @import ape
#' @author Ai-bing ZHANG,PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Zhang et al. 2015.Initinally created on 2014/8/2 11:23:41
#' @description kohonenClass2cclust transfers  "kohonen" class into "cclust" class in order to calculate cclustIndex using the function cclustIndex() which only takes an object of class "cclust" class as input.
#' @details only related data to calculating cclustIndex() were transferred!!!
#' @examples
#' require(kohonen)
#' data(wines)
#' set.seed(7)
#' training <- sample(nrow(wines), 120)
#' Xtraining <- scale(wines[training, ])
#' Xtest <- scale(wines[-training, ],
#' center = attr(Xtraining, "scaled:center"),
#'  scale = attr(Xtraining, "scaled:scale"))
#'  x <- som(Xtraining, grid = somgrid(5, 5, "hexagonal"))
#'  cl<-kohonenClass2cclust(x)
#'  class(cl)





library(testthat)
context("kohonenClass2cclust: transfers \"kohonen\" class into \"cclust\" class")
test_that("transfers \"kohonen\" class into \"cclust\" class",{
  require(kohonen)
  data(wines)
  set.seed(7)
  training <- sample(nrow(wines), 120)
  Xtraining <- scale(wines[training, ])
  Xtest <- scale(wines[-training, ],
                 center = attr(Xtraining, "scaled:center"),
                 scale = attr(Xtraining, "scaled:scale"))
  
  x <- som(Xtraining, grid = somgrid(5, 5, "hexagonal"))
  cl<-kohonenClass2cclust(x)
  
  expect_that(class(cl),
              equals("cclust")) 
  
})





