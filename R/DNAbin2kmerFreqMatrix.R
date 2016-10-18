#' Calculation of Kmer Frequency Matrix from DNAbin for Both Reference and Query Sequences
#' 
#' @description Calculation of kmer frequency matrices from DNAbin for both reference and query sequences.
#' 
#' @param  ref Object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  que Object of class "DNAbin", which needs to be inferred.
#' @param  kmer a numeric to indicate the length of kmer used.
#' @return kmer frequency matrices for both ref and que sequences, but only based on kmers found in ref!!!
#' new kmers in que will be ignored.
#' @keywords DNAbin2kmerFreqMatrix
#' @export 
#' 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA. zhangab2008(at)mail.cnu.edu.cn
#' @references 
#' zhangab2008(at)mail.cnu.edu.cn
#' 
#' 
#' 
#' 
#' 
#' @examples
#' data(TibetanMoth) 
#' ref<-as.DNAbin(as.character(TibetanMoth[1:50,]))
#' que<-as.DNAbin(as.character(TibetanMoth[51:60,]))
#' out<-DNAbin2kmerFreqMatrix(ref,que,kmer=3)
#' out
#' 
#' 

DNAbin2kmerFreqMatrix<-function(ref,que,kmer=kmer){ 
  ### return kmer frequency matrices for both ref and que sequences, but only based on kmers found in ref!!!
  ### new kmers in que will be ignored
  #require(seqinr)
  
  NAMES<-function(seqs){
    
    if(mode(seqs)=="raw"){
      SeqNames<-attr(seqs,"dimnames")[[1]]
    }else{ ### mode(seqs)=="list"
      SeqNames<-names(seqs)  
    }
    
    #names(SeqNames)<-NULL
    
    if(length(SeqNames)==0) {
      stop("the mode(seqs) is wrong!")
    }else{
      return(SeqNames) }
  }
  
  ##########
  c2s<-function (chars = c("m", "e", "r", "g", "e", "d")) {
    return(paste(chars, collapse = ""))
  }
  
  ##########
  
  
  
  
  ### 1. check the format of input arguments 
  if (class(ref)!="DNAbin")
    stop("seqs should be in DNAbin format!")
  
  if (class(kmer)!="integer")  kmer<-as.integer(kmer)
  
  ### 2. seek unique.kmer.vector for all seqs: u.s
  seqs.as.char<-as.character(ref)
  seqs.as.char2<-as.character(que)
  
  ifelse(is.vector(seqs.as.char),seqs.as.str.vector<-lapply(seqs.as.char,FUN=c2s),seqs.as.str.vector<-apply(seqs.as.char, MARGIN=1,FUN=c2s))
  ifelse(is.vector(seqs.as.char2),seqs.as.str.vector2<-lapply(seqs.as.char2,FUN=c2s),seqs.as.str.vector2<-apply(seqs.as.char2, MARGIN=1,FUN=c2s))
  
  
  seqs.unique.as.str.vector<-unique(seqs.as.str.vector)
  s<- seqs.unique.as.str.vector
  
  u.s0<-unique(substring(s, 1, kmer)) ### check along the column first!
  u.s<-u.s0 
  

  for (i in 2:(max(nchar(s))-kmer+1)){ 
    # i<-2
    #u.s<-u.s0   ### unique str
    u.s<-c(u.s,unique(substring(s, i, kmer+i-1)))
    u.s<-unique(u.s)
    n.char<-nchar(u.s)
    size.kmer.exact<-n.char==kmer ### logic to remove short kmer!
    u.s<-subset(u.s,size.kmer.exact)
    
    
  } ### the end of for-loop
  
  ### do some cleaning by removing IUPAC codes, "-"
  
  mpattern<-"-+[a-z]*"
  u.s<-gsub(mpattern,NA,u.s) # strings with "-"
  mpattern<-"[rymkswhbvdn]+"
  u.s<-gsub(mpattern,NA,u.s) # strings with "-"
  u.s<- u.s[!is.na(u.s)]
    ### 3.  calculate kmer frequency for each sequence in ref 
  #b<-gregexpr(u.s[1],seqs.as.str.vector[1])
  kmer.freq.one.seq<-sapply(u.s,gregexpr,seqs.as.str.vector[1])
  kmer.freq.one.seq2<-numeric(length(kmer.freq.one.seq))
  for (i in 1:length(kmer.freq.one.seq)){
    #kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]])
    ifelse(kmer.freq.one.seq[[i]]==-1,kmer.freq.one.seq2[i]<-0,kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]]))
  }        
  kmer.freq.one.seq3<- kmer.freq.one.seq2
  
  for (k in 2:length(seqs.as.str.vector)){
    # k<-2
    kmer.freq.one.seq<-sapply(u.s,gregexpr,seqs.as.str.vector[k])
    #b<-lapply(u.s,gregexpr,seqs.as.str.vector[1])
    kmer.freq.one.seq2<-numeric(length(kmer.freq.one.seq))
    
    for (i in 1:length(kmer.freq.one.seq)){
      
      ifelse(kmer.freq.one.seq[[i]]==-1,kmer.freq.one.seq2[i]<-0,kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]]))
    }
    
    kmer.freq.one.seq3<- c(kmer.freq.one.seq3,kmer.freq.one.seq2)
    
    
  } ### end of k-for-loop
  
  kmer.freq.matrix<-t(array(kmer.freq.one.seq3,dim=c(length(kmer.freq.one.seq),length(kmer.freq.one.seq3)%/%length(kmer.freq.one.seq))))
  #kmer.freq.matrix<-t(array(kmer.freq.one.seq3,dim=c(length(kmer.freq.one.seq3)%/%length(kmer.freq.one.seq),length(kmer.freq.one.seq)))) ### error!?
  #rownames(kmer.freq.matrix)<-rownames(ref)
  rownames(kmer.freq.matrix)<-NAMES(ref)
  ### 4.  calculate kmer frequency for each sequence in que 
  kmer.freq.one.seq<-sapply(u.s,gregexpr,seqs.as.str.vector2[1])
  
  kmer.freq.one.seq2<-numeric(length(kmer.freq.one.seq))
  for (i in 1:length(kmer.freq.one.seq)){
    #kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]])
    ifelse(kmer.freq.one.seq[[i]]==-1,kmer.freq.one.seq2[i]<-0,kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]]))
  }        
  kmer.freq.one.seq3<- kmer.freq.one.seq2
  for (k in 2:length(seqs.as.str.vector2)){
    # k<-2
    kmer.freq.one.seq<-sapply(u.s,gregexpr,seqs.as.str.vector2[k])
    #b<-lapply(u.s,gregexpr,seqs.as.str.vector[1])
    kmer.freq.one.seq2<-numeric(length(kmer.freq.one.seq))
    for (i in 1:length(kmer.freq.one.seq)){
      ifelse(kmer.freq.one.seq[[i]]==-1,kmer.freq.one.seq2[i]<-0,kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]]))
    }
    
    kmer.freq.one.seq3<- c(kmer.freq.one.seq3,kmer.freq.one.seq2)
  } ### end of k-for-loop
  
  kmer.freq.matrix2<-t(array(kmer.freq.one.seq3,dim=c(length(kmer.freq.one.seq),length(kmer.freq.one.seq3)%/%length(kmer.freq.one.seq))))
  #kmer.freq.matrix<-t(array(kmer.freq.one.seq3,dim=c(length(kmer.freq.one.seq3)%/%length(kmer.freq.one.seq),length(kmer.freq.one.seq)))) ### error!?
  
  #rownames(kmer.freq.matrix2)<-rownames(que)
  rownames(kmer.freq.matrix2)<-NAMES(que)
  
  out<-list(unique.str=u.s,kmer.Freq.ref=kmer.freq.matrix, kmer.Freq.que=kmer.freq.matrix2)
  #cat(kmer.freq.matrix)
  #return(kmer.freq.matrix)
  return(out)
}










