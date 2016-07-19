#' Fuzzy Membership Function Value
#' 
#' @description Calculation fuzzy membership function value given a distance from query to a potenial
#' species, maximual intraspecific variation of the potential species theta1, and minimal interspecific
#' distance (here, the distance between the potential species and its nearest neighbor theta2) 
#' (fuzzy-set based method, Zhang et al. 2012), different definition of distances could also be used.
#' 
#' @param  xtheta12 a numerical vector containing three elements, a distance from query to a potenial species,
#' maximual or sd of intraspecific variation of the potential species theta1,minimal or mean interspecific
#' distance.
#' @return a numeric between 0 and 1.
#' @keywords FMF
#' @export 
#' @author Ai-bing ZHANG, Zhi-yong SHI. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @note different definitions of distances could also be used.
#' @references 
#' Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach 
#' to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
#' @examples
#' 
#' xtheta12<-c(0.6289163,0.1465522,0.6379375) 
#' FMF.out<-FMF(xtheta12)
#' FMF.out

FMF<-function(xtheta12){
  ###  
  xtheta12<-as.numeric(xtheta12)
  
  if (class(xtheta12)!="numeric" ||length(xtheta12)!=3) 
    stop("input should be a numeric vector with length of 3!!!")
  
  x<-xtheta12[1]
  theta1<-xtheta12[2]
  theta2<-xtheta12[3]
  
  if (x<=theta1) FMF<-1
  if (x>theta1 && x<=(theta2+theta1)/2) FMF<-1-2*((x-theta1)/(theta2-theta1))^2
  if (x<=theta2 && x>(theta2+theta1)/2) FMF<-2*((x-theta2)/(theta2-theta1))^2
  if (x>=theta2) FMF<-0
  return(FMF)
}












