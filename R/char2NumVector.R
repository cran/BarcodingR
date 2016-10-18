#' Character to Integer Vector
#' 
#' @description Conversion from a character vector to an integer vector.
#' 
#' @param  c character vector.
#' @return an integer vector.
#' @keywords char2NumVector
#' @export 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA. zhangab2008 (at) mail. cnu. edu.cn.
#' @references 
#' zhangab2008 (at) mail. cnu. edu. cn.
#' 
#' @examples
#' c<-c("a","a","b")
#' num<-char2NumVector(c)
#' num
#' 
char2NumVector<-function(c){
  
  if (class(c)!="character") c<-as.character(c)
  
  #stop("invalid input format! input should be character vector!")
  c<-as.factor(c)
  
  level.c<-levels(c)
  #b<-dim(ecol.sample.list)[1]
  levels(c)<-seq(1:length(level.c))
  
  c<-as.numeric(c)
  
  return(c)}


