#' Consensus Identification
#' 
#' @description Make consensus for identifications from two or more methods, usually for 
#' a set of query sequences.
#' @param   identifiedBy2orMore an object of class "data.frame", containing (queIDs, as rownames), identifiedByMethod1,identifiedByMethod2,and so on. 

#' @return a data frame with consensus.identification, and corresponding votes.
#' @keywords consensus.identify
#' @export 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @note Suitable for case where a set of queries were identified by more than two methods.
#' @examples
#' 
#' queIDs<-c("q1","q2","q3")
#' 
#' bp<-c("sp1","sp1","sp1")
#' bpk<-c("sp1","sp1","sp2")
#' bayes<-c("sp2","sp1","sp3")
#' fuzzyID<-c("sp1","sp1","sp2")
#' identifiedBy2orMore<-data.frame(bp=bp,bpk=bpk,bayes=bayes,fuzzyID=fuzzyID)
#' rownames(identifiedBy2orMore)<-queIDs<-c("q1","q2","q3")
#' ccs<-consensus.identify(identifiedBy2orMore)



consensus.identify<-function (identifiedBy2orMore){
  
  concensus<-function(ids){ ### a character string indicating different identification for
                        ### by different methods
    ids<-as.factor(ids)
    
    t.ids<-table(ids)
     
    ifelse(max(t.ids)!=1,
           ccs<-levels(ids)[as.numeric(which.max(t.ids))],
           ccs<-"ambigous identification"
          )
    
    names(ccs)<-max(t.ids)
    return(ccs)
  } ### the end of concensus
  concensus2<-function(ids){ ### a character string indicating different identification for
    ### by different methods
    ids<-as.factor(ids)
    
    t.ids<-table(ids)
    
    #ifelse(max(t.ids)!=1,
    #       ccs<-levels(ids)[as.numeric(which.max(t.ids))],
    #       ccs<-"ambigous identification"
   # )
    
    #names(ccs)<-max(t.ids)
    return(max(t.ids))
  } 
 
  ccs<-apply(identifiedBy2orMore,MARGIN=1,FUN=concensus)
  ccs2<-apply(identifiedBy2orMore,MARGIN=1,FUN=concensus2)
  
  #ccs<-data.frame(queIDs=rownames(identifiedBy2orMore),concensus.id=ccs,votes=ccs2)
  ccs<-data.frame(concensus.id=ccs,votes=ccs2)
  

  return(ccs)
  

} 
### the end of the function

### how to call the function:

### ccs<-concensus.identify(identifiedBy2orMore)
### ccs
