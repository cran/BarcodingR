#' Optimize kmer Length
#' 
#' @description Optimize kmer length by trying kmers which length is in the range from 1 to max.kmer. 
#' The optimal kmer will have maximal species identification success rate.
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  max.kmer a numeric to indicate the length of maximal kmer.
#' 
#' @return a numeric indicating the optimal kmer in the range examined.
#' @keywords optimize.kmer
#' @export 
#' @author Ai-bing ZHANG, Cai-qing YANG, Meng-di HAO, CNU, Beijing, CHINA.
#' @references zhangab2008 (at) mail. cnu. edu. cn/zhangab2008 (at) gmail.com.
#' @examples
#' 
#' data(TibetanMoth) 
#' ref<-TibetanMoth[1:10,]
#' optimial.kmer<-optimize.kmer(ref,max.kmer=5)
#' 
#'
#'
optimize.kmer<-function (ref,max.kmer=max.kmer){
  #require(ape)
  set.seed(7)
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
  }### the end of the function
  ref.IDs<-NAMES(ref)
  #ref<-del.gaps(ref)
  #names(ref)<-ref.IDs
  
  morph.spe<-gsub(".+,","",ref.IDs) # remove sequence ID before ","  
  Spp2<-as.factor(morph.spe)
  
  strings.equal<-function(str1,str2){ifelse(str1==str2,1,0)}
  DNAbin2kmerFreqMatrix<-function(ref,kmer=kmer){ 
    ### return kmer frequency matricies for both ref and que sequences, but only based on kmers found in ref!!!
    ### new kmers in que will be ignored
    #require(seqinr)
    
    ##########
    c2s<-function (chars = c("m", "e", "r", "g", "e", "d")) {
      return(paste(chars, collapse = ""))
    }
    
    ##########
    
    #kmer<-1
    ### 1. check the format of input arguments 
    if (class(ref)!="DNAbin")
      stop("seqs should be in DNAbin format!")
    
    if (class(kmer)!="integer")  kmer<-as.integer(kmer)
    
    ### 2. seek unique.kmer.vector for all seqs: u.s
    seqs.as.char<-as.character(ref)
    #seqs.as.char2<-as.character(que)
    
    ifelse(is.vector(seqs.as.char),seqs.as.str.vector<-lapply(seqs.as.char,FUN=c2s),seqs.as.str.vector<-apply(seqs.as.char, MARGIN=1,FUN=c2s))
    #ifelse(is.vector(seqs.as.char2),seqs.as.str.vector2<-lapply(seqs.as.char2,FUN=c2s),seqs.as.str.vector2<-apply(seqs.as.char2, MARGIN=1,FUN=c2s))
    
    
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
    
    out<-list(unique.str=u.s,kmer.Freq.ref=kmer.freq.matrix)
    
    #cat(kmer.freq.matrix)
    return(kmer.freq.matrix)
    
    #return(out)
    
  }#
  
  success.rates<-numeric(max.kmer)
  
  for (i in 1:max.kmer){
    
    # i<-1
    
    kmer.freq.ref<-DNAbin2kmerFreqMatrix(ref,kmer=i)
    
    
    knn1<-knn(kmer.freq.ref, kmer.freq.ref, cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    
    #knn1<-knn1(kmer.freq.ref, kmer.freq.ref, cl=Spp2)
    #attributes(.Last.value)
    #attributes(knn1)
    
    #knn1
    spe.morph<-as.character(Spp2)
    spe.Identified<-as.character(knn1)
    
    spe.morph.Identified<-data.frame(spe.morph,spe.Identified,stringsAsFactors=TRUE)
    
    matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])
    
    matches<-colSums(matches,dims = 1)
    
    success.rates[i]<-matches[1]/matches[2]
    
    #cat("i:",i,"\n")
    #cat("spe.morph:",spe.morph,"\n")
    #cat("spe.Identified:",spe.Identified,"\n")
    
    #cat("success.rates[i]:",success.rates[i],"\n")
  }
  
  kmer.best<-which.max(success.rates)
  
  ### plot start
  plot(1:max.kmer,success.rates,type="h",main="Success rates of spe identification with different length of kmer (ref)",
       xlab="k (length of kmer)", ylab="Success rates of spe identification")
  axis(1,kmer.best, paste("optimum",kmer.best,sep="\n"),col="red",font=2,col.axis="red")
  points(kmer.best,max(success.rates),pch=16,col="red",cex=1.5)
  ### plot end
  
 
  success.rates.ref<-max(success.rates)
  names(success.rates.ref)<-NULL
  out<-c(kmer.best,success.rates.ref)
  
  
  
  return(out)
  #return(kmer.best)
  
}
























