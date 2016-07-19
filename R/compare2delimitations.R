#' Comparision between two Delimitations
#' 
#' @description Comparision between two delimitations of a group of samples, for instance,
#' traditionally morphological delimitation and molecular delimitation (MOTU).
#' 
#' @param  deli1 a character array (vector),containing a set of, for example, morphological identification (species names), to compare with
#' @param  deli2 a character array (vector),containing a set of, molecular delimitation (MOTU). 
#' @return a list containing the adjusted Rand index comparing the two partitions (a scalar). This index has zero expected value in the case
#'  of random partition, and it is bounded above by 1 in the case of perfect agreement between two partitions; the numbers of matches,
#'  splits,merges, and corresponding percentage.
#' @keywords compare2delimitations
#' @export 
#' 
#' @note This is for the same set of samples with two partitions/delimitations.
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA.
#' @references
#' L. Hubert and P. Arabie (1985) Comparing Partitions, Journal of the Classification 2:193-218.
#' @examples
#' 
#' 
#' 
#' 
#' 
#' deli1<-c(1,1,1,1,1,1)
#' deli2<-c(1,1,2,1,1,3)
#' out<-compare2delimitations(deli1,deli2)
#' out


compare2delimitations<-function(deli1,deli2){
  char2NumVector<-function(c){
    
    if (class(c)!="character") c<-as.character(c)
    
    #stop("invalid input format! input should be character vector!")
    c<-as.factor(c)
    
    level.c<-levels(c)
    #b<-dim(ecol.sample.list)[1]
    levels(c)<-seq(1:length(level.c))
    
    c<-as.numeric(c)
    
    return(c)}
  
  ##########
  adjustedRandIndex<-function (x, y){
    x <- as.vector(x)
    y <- as.vector(y)
    if (length(x) != length(y)) 
      stop("arguments must be vectors of the same length")
    tab <- table(x, y)
    if (all(dim(tab) == c(1, 1))) 
      return(1)
    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2)) - a
    c <- sum(choose(colSums(tab), 2)) - a
    d <- choose(sum(tab), 2) - a - b - c
    ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                       a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    return(ARI)
  }
  
  ##########
  
  
  
  #U<-deli1  ### morphological classification, being compared to
  #V<-deli2
  
  if(length(deli1)!=length(deli2))
    stop("deli1 and deli2 shoud be same in length!")
  
  twoPartitions<-cbind(deli1,deli2)  
  no<-1:length(deli1)
  
  rownames(twoPartitions)<-paste("sample",no,sep="")
  
  U<-twoPartitions[,1]  ### morphological classification, being compared to
  V<-twoPartitions[,2]
  
  
  
  no.partion.U<-length(levels(factor(U)))
  no.partion.V<-length(levels(factor(V)))  
  
  ### to generate contingency table of U and V
  ### ContingencyTableUV:ContableUV
  k=1 ### the total number of elements in ContableUV
  #rm (no.in.common)
  #no.in.common <- c(0) ### initilize the vector!??
  no.in.common <- NULL ### initilize the vector!
  
  for (j in 1: no.partion.V){
    
    for (i in 1:no.partion.U){
      
      sub.U<-U[U==levels(factor(U))[i]]
      
      names.sub.U<-names(sub.U)
      
      sub.V<-V[V==levels(factor(V))[j]]
      
      names.sub.V<-names(sub.V)
      
      
      no.in.common[k]<-length(intersect(names.sub.U,names.sub.V))
      #ContableUV[i,j]=length(intersect(names.sub.U,names.sub.V))
      #no.in.common[k]
      k=k+1
    }
  }
  
  ContableUV<-array(no.in.common,c(no.partion.U,no.partion.V))
  #ContableUV<-array(no.in.common,c(no.partion.V,no.partion.U))
  #ContableUV
  
  ### to count the number and frequency of MATCH, SPLIT, MERGE, MIXTURE(?) of V (e.g. MOTUS)
  ###         referred to U (morpho-partition)
  MATCH<-0
  SPLIT<-0
  indexTrue<-0
  for (j in 1:no.partion.U) {
    #j<-3
    cat("j=",j)
    cat("\n")
    #x<-ContableUV[,j]
    x<-ContableUV[j,]
    if (length(x[x>0])==1&&length(x)!=1) {
      a<-x>0
      for (k in 1:length(a)) {
        #k<-3
        #ifelse(a[k]==TRUE,indexTrue<-k,break)
        if(a[k]==TRUE) indexTrue<-k
      }
      x2<-ContableUV[,indexTrue]
      
      if (length(x2[x2>0])==1) MATCH<-MATCH+1
      
      #else SPLIT<-SPLIT+1
      
    }
    else SPLIT<-SPLIT+1
    
    
  }
  #freq.MATCH=MATCH/no.partion.V
  #freq.SPLIT=SPLIT/no.partion.V
  freq.MATCH=MATCH/no.partion.U
  freq.SPLIT=SPLIT/no.partion.U
  
  
  ### to count MERGE
  MERGE<-0
  for (i in 1:no.partion.V) {
    #x<-ContableUV[i,]
    x<-ContableUV[,i]
    #if (length(x[x>0])==1) MATCH<-MATCH else MERGE<-MERGE+1
    if (length(x[x>0])>1) MERGE<-MERGE+1
  }
  
  freq.MERGE=MERGE/no.partion.V
  
  # freq.mixture<-1-freq.MATCH-freq.SPLIT-freq.MERGE
  
  ### to calculate ARI:The adjusted Rand index comparing the two partitions.
  deli1<-char2NumVector(deli1)
  deli2<-char2NumVector(deli2)
  
  ARI<-adjustedRandIndex(deli1, deli2)
  
  
  #OUT<-list(Match=MATCH,Split=SPLIT,Merge=MERGE,Freq.match=freq.MATCH,Freq.split=freq.SPLIT,Freq.merge=freq.MERGE,Freq.mixture=freq.mixture,ContableUV=ContableUV)
  OUT<-list(adjustedRandIndex=ARI,Match=MATCH,Split=SPLIT,Merge=MERGE,Freq.match=freq.MATCH,Freq.split=freq.SPLIT,Freq.merge=freq.MERGE,ContableUV=ContableUV)
  
  return(OUT)
  
  
  
}


















