#' Bp Barcoding Species Identify using Kmer
#' 
#' @description Species identification using BP-based method for both protein-coding
#' barcodes, for instance, COI, and non-coding barcodes, such as, ITS, using kmer statistics.
#' 
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  que object of class "DNAbin", which needs to be inferred.
#' @param  kmer a numeric indicating the length of kmer used.
#' @param  UseBuiltModel logic value to indicate whether a built model is used 
#' or not.
#' @param  lr parameter for weight decay. Default 5e-5.
#' @param  maxit maximum number of iterations. Default 1e+6.
#' @return a list containing model parameters used, species identification success rates using references,
#' query sequences, species inferred, and corresponding confidence levels (bp probability for BP-based method).
#' @keywords bbsik
#' @export 
#' @import nnet

#' @author Ai-bing ZHANG, Meng-di HAO, Cai-qing YANG, CNU, Beijing, CHINA. zhangab2008 (at) mail. cnu. edu.cn
#' @references 
#' Zhang, A. B., D. S. Sikes, C. Muster, S. Q. Li. (2008). Inferring Species Membership 
#' using DNA sequences with Back-propagation Neural Networks. Systematic Biology, 57(2):202-215. 
#' 
#' 
#' @examples
#' data(TibetanMoth) 
#' ref<-as.DNAbin(as.character(TibetanMoth[1:50,]))
#' que<-as.DNAbin(as.character(TibetanMoth[51:60,]))
#' out<-bbsik(ref, que, kmer = 1, UseBuiltModel = FALSE)
#' out
#' out$convergence
#' out$success.rates.ref
#' 
#' data(pineMothITS2) 
#' ref<-pineMothITS2
#' que<-pineMothITS2
#' out<-bbsik(ref, que, kmer = 1, UseBuiltModel = FALSE)
#' out
#' out$convergence
#' out$success.rates.ref
#' 
#' 
#' 

bbsik<-function (ref, que, kmer = kmer, UseBuiltModel = FALSE,lr=5e-5, maxit=1e+6) {
  
  ### functions used!
  DNAbin2kmerFreqMatrix2<-function(ref,que,kmer=kmer){ 
    ### return kmer frequency matrices for both ref and que sequences, but only based on kmers found in ref!!!
    ### new kmers in que will be ignored
    #require(seqinr)
    
    #kmer<-1
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
    
    #length(u.s)
    
    
    
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
    
    rownames(kmer.freq.matrix)<-rownames(ref)
    
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
    
    rownames(kmer.freq.matrix2)<-rownames(que)
    out<-list(unique.str=u.s,kmer.Freq.ref=kmer.freq.matrix, kmer.Freq.que=kmer.freq.matrix2)
    
    #cat(kmer.freq.matrix)
    #return(kmer.freq.matrix)
    
    return(out)
    
  }
  
  char2NumVector<-function(c){
    
    if (class(c)!="character") c<-as.character(c)
    
    #stop("invalid input format! input should be character vector!")
    c<-as.factor(c)
    
    level.c<-levels(c)
    #b<-dim(ecol.sample.list)[1]
    levels(c)<-seq(1:length(level.c))
    
    c<-as.numeric(c)
    
    return(c)}
  
  strings.equal<-function(str1,str2){ifelse(str1==str2,1,0)}
  ##########
  c2s<-function (chars = c("m", "e", "r", "g", "e", "d")) {
    return(paste(chars, collapse = ""))
  }
  ##########
  
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
  
  

  if(UseBuiltModel == FALSE){
  
    
    sampleSpeNames<-attr(ref,"dimnames")[[1]]
    if(length(sampleSpeNames)==0) sampleSpeNames<-names(ref)
    
    
    queIDs<-attr(que,"dimnames")[[1]]
    if(length(queIDs)==0) queIDs<-names(que)
    
    
    mpattern<-".+,"
    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    #Spp
    #length(unique(Spp))
    sampleSpeNames2<-class.ind(Spp)
    
    #if (kmer>0.05*len.shortest.seq)
    #  stop("kmer is too large, it will take lots of time to run! a vaule less than 10 is suggested!")
    out1<-DNAbin2kmerFreqMatrix2(ref,que,kmer=kmer)
    
    
    unique.str.ref<-out1$unique.str
    
   # bp<-function(ref1,sampleSpeNames2,que1,unique.str.ref,lr=5e-5, maxit=1e+6){
      bp<-function(ref1,sampleSpeNames2,que1,lr=5e-5, maxit=1e+6){
      #unique.str.ref<-out1$unique.str
      
      
      ### standard
      center.ref1<-apply(ref1,MARGIN=2,FUN=mean)
      sd.ref1<-apply(ref1,MARGIN=2,FUN=sd)
      ref1<-scale(ref1,scale=F)
      que1<-scale(que1,center = center.ref1, scale=F)
      #que1<-scale(que1,center = center.ref1, scale=sd.ref1)
      
      ref1_tmp<-ref1[!is.na(ref1)]
      range1<-1./max(abs(ref1_tmp))
      # range1

      n.hidden<-ceiling(log2(dim(ref1)[1]))
         
      nnet.trained <- nnet(ref1, sampleSpeNames2, size = n.hidden, rang = range1,### 0.5
                           entropy = FALSE,
                           decay = lr, maxit = maxit,#decay = 5e-5, maxit = 1e+6,
                           abstol = 1.0e-8, 
                           reltol = 1.0e-8,
                           MaxNWts = 20000) ### decay = 5e-5
      
      spe.inferred0<-predict(nnet.trained, que1)
      spe.inferred<-spe.inferred0
      #spe.inferred[spe.inferred>=0.95]<-1
      #spe.inferred[spe.inferred<0.95]<-0
     
      #colnames(spe.inferred)
      #which.max(spe.inferred[1,])
      inferred<-apply(spe.inferred,1,FUN=which.max)
      inferred.prob<-apply(spe.inferred,1,FUN=max)
      # length(inferred)
      
      #inferred
      #colnames(spe.inferred)[inferred]
      
      
      ############################################################
      ####### calculate model success rate using ref  start...
      ###########################################################
      spe.inferred0.ref<-predict(nnet.trained, ref1)
      spe.inferred.ref<-spe.inferred0.ref
      
      inferred.ref<-apply(spe.inferred.ref,1,FUN=which.max)
      inferred.prob.ref<-apply(spe.inferred.ref,1,FUN=max)
      # length(inferred)
      
      #inferred
      p.ref<-colnames(spe.inferred.ref)[inferred.ref]
      
      #p.ref[1]<-"abz"
      
      spe.morph.Identified<-data.frame(Spp,p.ref)
      
      matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])
      
      
      matches<-colSums(matches,dims = 1)
      
      success.rates.ref<-matches[1]/matches[2]
      names(success.rates.ref)<-NULL
      
      ############################################################
      ####### calculate model success rate using ref the end.
      ###########################################################
      
      
      output.identified<-data.frame(queIDs=queIDs,
                                    spe.Identified=colnames(spe.inferred)[inferred],
                                    bp.prob=inferred.prob)
      #output.identified
      
      rownames(output.identified)<-NULL
      
      
      
      
      
      out.bp<-list(summary.model=summary(nnet.trained),###convergence=nnet.trained$convergence,
                   convergence=nnet.trained$convergence,
                   success.rates.ref=success.rates.ref,
                   kme.ref=kmer,
                   output_identified=output.identified,
                   nnet.trained=nnet.trained,
                   center.ref1=center.ref1
                   )
      
      
      
      #current.wd<-getwd()
      #setwd
      
      
      
      return(out.bp)
      
      #return(output.identified)
      
    }
    
    output2.identified<-bp(ref1=out1$kmer.Freq.ref,
             sampleSpeNames2,
             que1=out1$kmer.Freq.que,
             #unique.str.ref=out1$unique.str,
             lr=5e-5, 
             maxit=1e+6)
    
    
    fileName<-"bbsik_tmp"
    #fileName<-paste("simulation",i,sep = "")
    fileName<-paste(fileName,".RData",sep = "")
    save(output2.identified,
      #nnet.trained,
      #   out.bp,
       #  success.rates.ref,
      #   kmer,
         unique.str.ref,
       #  center.ref1,
         file = fileName) ## #save(x, y, file = "xy.RData") 
    
      
   
  }else{
    
    if(!file.exists("bbsik_tmp.RData"))
    stop("file bbsik_tmp.RData not found in current directory!
         you need to move the file into current directory or rebuild the model!
         ")
    
    ### 4.  calculate kmer frequency for each sequence in que 
    load("bbsik_tmp.RData")
    center.ref1<-output2.identified$center.ref1
    nnet.trained<-output2.identified$nnet.trained
    success.rates.ref<-output2.identified$success.rates.ref
    queIDs<-attr(que,"dimnames")[[1]]
    if(length(queIDs)==0) queIDs<-names(que)
    #u.s<-unique.str.ref
    u.s<-unique.str.ref
    seqs.as.char2<-as.character(que)
    ifelse(is.vector(seqs.as.char2),seqs.as.str.vector2<-lapply(seqs.as.char2,FUN=c2s),seqs.as.str.vector2<-apply(seqs.as.char2, MARGIN=1,FUN=c2s))
    
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
      kmer.freq.one.seq2<-numeric(length(kmer.freq.one.seq))
      for (i in 1:length(kmer.freq.one.seq)){  
        ifelse(kmer.freq.one.seq[[i]]==-1,kmer.freq.one.seq2[i]<-0,kmer.freq.one.seq2[i]<-length(kmer.freq.one.seq[[i]]))
      }
      
      kmer.freq.one.seq3<- c(kmer.freq.one.seq3,kmer.freq.one.seq2)
    } ### end of k-for-loop
    
    kmer.freq.matrix2<-t(array(kmer.freq.one.seq3,dim=c(length(kmer.freq.one.seq),length(kmer.freq.one.seq3)%/%length(kmer.freq.one.seq))))
    rownames(kmer.freq.matrix2)<-rownames(que)
    if(length(rownames(kmer.freq.matrix2))==0) rownames(kmer.freq.matrix2)<-names(que)
      que1<-scale(kmer.freq.matrix2,center = center.ref1, scale=F)
      spe.inferred0<-predict(nnet.trained, que1)
      spe.inferred<-spe.inferred0
      inferred<-apply(spe.inferred,1,FUN=which.max)
      inferred.prob<-apply(spe.inferred,1,FUN=max)
      output.identified<-data.frame(queIDs=queIDs,
                                  spe.Identified=colnames(spe.inferred)[inferred],
                                  bp.prob=inferred.prob) 
    rownames(output.identified)<-NULL
    
    output2.identified<-list(summary.model=summary(nnet.trained),
                             convergence=nnet.trained$convergence,
                             success.rates.ref=success.rates.ref,
                             kme.ref=kmer,
                             output_identified=output.identified)
  }
  
  class(output2.identified) <- c("BarcodingR")
  
  return(output2.identified)
}  ### the end of the function!


