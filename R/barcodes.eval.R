#' Barcodes Evaluation
#' 
#' @description Evaluate two barcodes using species identification success rate critera.
#' 
#' @param  barcode1 object of class "DNAbin" based on barcode1, which contains taxon information.
#' @param  barcode2 object of class "DNAbin" based on barcode2, which contains taxon information.
#' @param  kmer1 a numeric to indicate the length of kmer1 for barcode1, the opitimal kmer could be found by the function
#' optimize.kmer() before running this function.
#' @param  kmer2 a numeric to indicate the length of kmer2 for barcode2, see above.
#' 
#' @return a list containing p_value of prop.test(), and so on.
#' @keywords barcodes.eval
#' @export 
#' @import ape
#' @import class
#' 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA.
#' @references zhangab2008(at)mail.cnu.edu.cn
#' @seealso prop.test()
#' @examples
#' 

#' data(TibetanMoth) 
#' barcode1<-as.DNAbin(as.character(TibetanMoth[1:30,]))
#' barcode2<-barcode1
#' b.eval<-barcodes.eval(barcode1,barcode2,kmer1=1,kmer2=3)
#' b.eval

#b.eval<-barcodes.eval(pineMothITS2,pineMothITS1,kmer1=3,kmer2=3)
#b.eval

#library(help=class)

barcodes.eval<-function (barcode1,barcode2,kmer1=kmer1,kmer2=kmer2){
  set.seed(7)
  ### general
  #barcode1.IDs<-rownames(barcode1)
  #barcode2.IDs<-rownames(barcode2)
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
  barcode1.IDs<-NAMES(barcode1)
  barcode2.IDs<-NAMES(barcode2)
  
  barcode1<-del.gaps(barcode1)
  barcode2<-del.gaps(barcode2)
  names(barcode1)<-barcode1.IDs
  names(barcode2)<-barcode2.IDs
  
  ### functions:
  strings.equal<-function(str1,str2){ifelse(str1==str2,1,0)}
  
  ##########
  c2s<-function (chars = c("m", "e", "r", "g", "e", "d")) {
    return(paste(chars, collapse = ""))
  }
  
  ##########
  
  
  DNAbin2kmerFreqMatrix<-function(ref,kmer=kmer){ 
    ### return kmer frequency matricies for both ref and que sequences, but only based on kmers found in ref!!!
    ### new kmers in que will be ignored
    #require(seqinr)
    
    #kmer<-3
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
    
    ifelse(length(names(ref))!=0,
           rownames(kmer.freq.matrix)<-names(ref),
           rownames(kmer.freq.matrix)<-rownames(ref))
    
    out<-list(unique.str=u.s,kmer.Freq.ref=kmer.freq.matrix)
    
    #cat(kmer.freq.matrix)
    return(kmer.freq.matrix)
    
    #return(out)
    
  }#
  

  ### marker1
  success.b1<-0
  morph.spe<-gsub(".+,","",barcode1.IDs) # remove sequence ID before ","
  no.samples.b1<-length(barcode1)
  barcode0<-barcode1
  
  for(i in 1:no.samples.b1){
    # i<-1
    barcode1<-barcode0
    
    #barcode1.minus1<-barcode1[-i,]
    #barcode1<-as.matrix(barcode1)
    #barcode1.minus1<-barcode1[-i,]
    morph.spe.minus1<-morph.spe[-i]
    #length(morph.spe.minus1)
    
    kmer.freq.b0<-DNAbin2kmerFreqMatrix(barcode1,kmer=kmer1)
    #kmer.freq.b1<-DNAbin2kmerFreqMatrix(barcode1.minus1,kmer=kmer1)
    test<-kmer.freq.b0[i,]
    kmer.freq.b0.minus1<-kmer.freq.b0[-i,]
    
    #cat("i:",i,"\n")
    #cat("size.kmer.freq.b1=",dim(kmer.freq.b0)[1],"\n")
    
    Spp2<-as.factor(morph.spe.minus1)
    
    knn1<-knn(kmer.freq.b0.minus1, 
              test, 
              cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    spe.Identified<-as.character(knn1)
    
    #success.b1<-strings.equal(as.character(Spp2)[i],spe.Identified)+ success.b1
    success.b1<-strings.equal(morph.spe[i],spe.Identified)+ success.b1
    
    
  }
    
  #success.b1
  
  
  
  ### marker2
  success.b2<-0
  morph.spe2<-gsub(".+,","",barcode2.IDs) # remove sequence ID before ","
  no.samples.b2<-length(barcode2)
  barcode0<-barcode2
  
  for(i in 1:no.samples.b2){
    # i<-1
    barcode1<-barcode0
    
    #barcode1.minus1<-barcode1[-i,]
    #barcode1<-as.matrix(barcode1)
    #barcode1.minus1<-barcode1[-i,]
    morph.spe.minus1<-morph.spe2[-i]
    #length(morph.spe.minus1)
    
    kmer.freq.b0<-DNAbin2kmerFreqMatrix(barcode1,kmer=kmer1)
    #kmer.freq.b1<-DNAbin2kmerFreqMatrix(barcode1.minus1,kmer=kmer1)
    test<-kmer.freq.b0[i,]
    kmer.freq.b0.minus1<-kmer.freq.b0[-i,]
    
    cat("i:",i,"\n")
    #cat("size.kmer.freq.b1=",dim(kmer.freq.b0)[1],"\n")
    
    Spp2<-as.factor(morph.spe.minus1)
    
    knn1<-knn(kmer.freq.b0.minus1, 
              test, 
              cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    spe.Identified<-as.character(knn1)
    
    #success.b1<-strings.equal(as.character(Spp2)[i],spe.Identified)+ success.b1
    success.b2<-strings.equal(morph.spe2[i],spe.Identified)+ success.b2
    
    
  }
  
  #success.b2
  
  ### prop.test()
  
  success<-c(success.b1,success.b2)
  total<-c(no.samples.b1,no.samples.b2)
  
  out<-prop.test(success,total,alternative ="greater")
  
  
  out2<-list(X_squared = out$statistic,
             p.value = out$p.value,
             estimate = out$estimate,
             conf.int = out$conf.int)
  
  
  return(out2)
  
  
}








