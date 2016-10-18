#' Species Identification Based on Fuzzy-set Method and kmer
#' @description Species identification based on fuzzy-set method (Zhang et al. 2012)and kmer.
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  que object of class "DNAbin", whose identities (species names) need to be inferred.
#' @param  kmer a numeric to indicate the length of maximum kmer to try in the range of 1 to kmer in case of
#' optimization = TRUE, otherwise, only a certain length of kmer is used.
#' @param optimization a character string, indicating whether different length of kmer (up to kmer) will be used or
#' just a specified length of kmer will be used.
#' @return a list indicating the identified species.
#' @keywords barcoding.spe.identify2
#' @export 
#' @import class
#' @import sp
#' @author Ai-bing ZHANG, Cai-qing YANG, Meng-di HAO, CNU, Beijing, CHINA, contact at zhangab2008 (at) mail. cnu. edu. cn.
#' @note read.dna() from package {ape} was used to obtain DNAbin object for unaligned non-coding barcodes. 
#'   
#' @references 
#' 
#' Zhang, A.B, Hao, M.D., Yang,C.Q., Shi, Z.Y. (2016). BarcodingR: an integrated R package for species identification using DNA barcodes. Methods in Ecology and Evolution. In press.
#' 
#' Jin,Q., H.L. Han, X.M. Hu, X.H. Li,C.D. Zhu,S. Y. W. Ho, R. D. Ward, A.B. Zhang . (2013).Quantifying Species Diversity with a DNA Barcoding-Based Method: Tibetan Moth Species (Noctuidae) on the Qinghai-Tibetan Plateau. PloS One 8: e644.
#' 
#' Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
#' 
#' Zhang, A. B., D. S. Sikes, C. Muster, S. Q. Li. (2008). Inferring Species Membership using DNA sequences with Back-propagation Neural Networks. Systematic Biology, 57(2):202-215. 
#' @examples
#' 
#' data(pineMothITS2) 
#' ref<-pineMothITS2
#' que<-ref
#' spe.id<-barcoding.spe.identify2(ref,que, kmer = 1, optimization = FALSE)
#' spe.id
#'  
barcoding.spe.identify2<-function (ref, que, kmer = kmer, optimization = TRUE) { # "fixedNumber","optimal"
  set.seed(7)
  ref.IDs<-rownames(ref)
  que.IDs<-rownames(que)
  
  if(length(ref.IDs)==0) ref.IDs<-names(ref)
  if(length(que.IDs)==0) que.IDs<-names(que)
  
  if(length(ref.IDs)!=0){
  que.IDs<-rownames(que)
  if(length(que.IDs)==0) que.IDs<-names(que)
  
  #ref<-del.gaps(ref)
  #que<-del.gaps(que)
  #names(ref)<-ref.IDs
  #names(que)<-que.IDs
  }
  
  ##########
  c2s<-function (chars = c("m", "e", "r", "g", "e", "d")) {
    return(paste(chars, collapse = ""))
  }
  
  ##########
  
  
  
  
  morph.spe<-gsub(".+,","",ref.IDs) # remove sequence ID before ","  
  Spp2<-as.factor(morph.spe)
  
  len.shortest.seq<-min(as.numeric(summary(ref)[,1]))
    
  #if (kmer>0.05*len.shortest.seq)
   # stop("kmer is too large, it will take lots of time to run! a vaule less than 10 is suggested!")
  
 # optimization <- TRUE
  #################################################
  ### functions used!
  #################################################
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
      
      if(length(rownames(kmer.freq.matrix))==0) rownames(kmer.freq.matrix)<-names(ref)
      
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
      if(length(rownames(kmer.freq.matrix2))==0) rownames(kmer.freq.matrix2)<-names(que)
      
      out<-list(unique.str=u.s,kmer.Freq.ref=kmer.freq.matrix, kmer.Freq.que=kmer.freq.matrix2)
      
      #cat(kmer.freq.matrix)
      #return(kmer.freq.matrix)
      
      return(out)
      
    }### the end of function
  
    strings.equal<-function(str1,str2){ifelse(str1==str2,1,0)}
  
    FMF<-function(xtheta12){
    ###  
    xtheta12<-as.numeric(xtheta12)
    
    if (class(xtheta12)!="numeric" ||length(xtheta12)!=3) 
      stop("input should be a numeric vector with length of 3!!!")
    
    x<-xtheta12[1]
    theta1<-xtheta12[2]
    theta2<-xtheta12[3]
    
    ##### test:
    #x<-0.6289163
    #theta1<-0.1465522
    #theta2<-0.6379375
    ##### test:
    
    
    
    if (x<=theta1) FMF<-1
    if (x>theta1 && x<=(theta2+theta1)/2) FMF<-1-2*((x-theta1)/(theta2-theta1))^2
    if (x<=theta2 && x>(theta2+theta1)/2) FMF<-2*((x-theta2)/(theta2-theta1))^2
    if (x>=theta2) FMF<-0
    return(FMF)
  }
    
    FMFtheta12<-function(Ref){
    
    # if (class(seqsRef)!="matrix")
    
    # stop("input should be an object of matrix!")
    
    
    ### 2.1 dealing with input DNA data! 
    ### 2.2 calculate species center vectors  
     # rm(ref)
     # Ref<-pineMothITS2
    morph.spe<-gsub(".+,","",rownames(Ref)) # remove sequence ID before ","   
    if(length(morph.spe)==0) morph.spe<-gsub(".+,","",names(Ref))
    #no.morph.spe<-length(unique(morph.spe))
    #species.centers<-aggregate(scale(Ref),by=list(morph.spe),FUN="mean")
    species.centers<-aggregate(Ref,by=list(morph.spe),FUN="mean")
    
    
    
    list.spe<-species.centers[,1]
    species.centers<-species.centers[,-1]
    
    ### 2.3  seek NN for PS (all species in this case!)
    ### 2.3.1.calculate pair distance of species centers  
    units.dist<-dist(species.centers, method = "euclidean", diag = F, upper = T, p = 2)
    units.dist0<-units.dist
    units.dist<-as.matrix(units.dist) ### important!
    #units.dist0<-units.dist
    
    for (i in 1: nrow(units.dist)) {units.dist[i,i]<-NA}
    
    
    #####
    ### 2.4. look for elements (their indices) with minimal distance to each other
    index1<-numeric(length(unique(morph.spe)))
    index2<-index1
    min.dist<-index1
    
    
    for (i in 1:nrow(units.dist)){
      # i<-1
      index1[i]<-i
      
      b<-which.min(units.dist[i,])
      
      if (length(b)==0) 
      {index2[i]<-NA
       min.dist[i]<-NA}
      else {index2[i]<-b
            min.dist[i]<-min(units.dist[i,],na.rm=T)}
    } ### for loop
    
    
    pairs<-rbind(index1,index2)
    pairs<-t(pairs)
    #class(pairs)
    
    pairs<-subset(pairs,subset=!is.na(pairs[,2]))
    
    
    
    theta1.tmp<-numeric(length(unique(morph.spe)))
    
    
    Spp<-morph.spe  ### seqsRef$unit.classif
    #seqs<-scale(digitized.locus) ### seqsRef$data
    #dim(seqs)
    uniSpeNames<-unique(Spp)
    
    
    #table(Spp)
    
    f<-factor(Spp)
    
    
    #####
    ###################################
    ### popSize calculation
    ###################################
    popSize.PS<-as.numeric(table(Spp))
    
    
    ###################################
    ### theta1,2 calculation
    ###################################
    
    ### for i<-1
    
    #source("eucl.dist.two.vect.R")
    
    seqInOneSpe<-Ref[grep(levels(f)[1], Spp, value = FALSE,fixed = TRUE),]
    #dim(seqInOneSpe)==NULL
    ifelse(popSize.PS[1]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean)) ###
    ifelse(popSize.PS[1]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) ###  centroid.spe0 - 
    
    length.sites<-length(centroid.spe)
    ifelse(popSize.PS[1]==1,intra.dist<-0,intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe)) ### to the centroid
    ifelse(popSize.PS[1]==1,theta1.tmp[1]<-0,theta1.tmp[1]<-sd(intra.dist)) ### theta1 is sligthly different from 
    #ifelse(popSize.PS[1]==1,theta1.tmp[1]<-0,theta1.tmp[1]<-max(intra.dist)) ### theta1 is sligthly different from 
    
    ######
    
    
    
    for(i in 2:length(levels(f))){
      #i<-2
      #   i<-9
      
      #i=3
      cat(paste("i=",i),"\n")
      # cat("\n")
      
      
      #seqInOneSpe<-sDNAbin[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),]
      seqInOneSpe<-Ref[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),] 
      ifelse(popSize.PS[i]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) 
      #ifelse(popSize.PS[i]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean))
      ifelse(popSize.PS[i]==1,centroid.spe<-c(centroid.spe,seqInOneSpe),centroid.spe<-c(centroid.spe,apply(seqInOneSpe,2,mean)))
      
      
      #length(seqInOneSpe) 
      
      #ifelse(dim(seqInOneSpe)==NULL,theta1.tmp[i]<-0,theta1.tmp[i]<-max(dist(seqInOneSpe))) # error!
      #ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-max(dist(seqInOneSpe)))
      
      ifelse(popSize.PS[i]==1,intra.dist<-eucl.dist.two.vect(seqInOneSpe,centroid.spe0),intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0))
      
      
      #intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0) ### to the centroid 
      #ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-max(intra.dist)) ### theta1 is sligthly different from 
      ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-sd(intra.dist)) ### theta1 is sligthly different from 
      
      
      #theta1.tmp[i]<-max(dist(seqInOneSpe))
      
      
    }
    
    
    
    
    centroid.spe.matrix<-t(array(centroid.spe,dim=c(length.sites,length(centroid.spe)%/%length.sites)))
    
    
    ### 1.2 calculate theta2 for each pairt of species
    
    
    ###################
    #codes<-out.somu$out.som.unique$codes  ### seqsRef$codes
    codes<-centroid.spe.matrix
    
    # dim(codes)
    
    
    theta12.all<-data.frame(list.spe=list.spe,PS=pairs[,1],NN=pairs[,2])
    
    #source("eucl.dist.two.vect.R")
    theta2.tmp<-numeric(dim(pairs)[1])
    
    
    for (i in 1:dim(pairs)[1]){
      # i<-1
      #v1<-codes[i, pairs[i,1]]
      #v2<-codes[i, pairs[i,2]]
      
      v1<-codes[pairs[i,1],]
      v2<-codes[pairs[i,2],]
      
      
      theta2.tmp[i]<-eucl.dist.two.vect(v1,v2)
    }
    
    
    theta12.all$theta1<-with(theta12.all,theta1.tmp)
    theta12.all$theta2<-with(theta12.all,theta2.tmp)
    
    theta12.all$popSize.PS<-with(theta12.all,popSize.PS)
    
    #theta12.all
    return(theta12.all) 
    
    #class(theta12.all)
    
  } ### the end of the function
    
    eucl.dist.two.vect<-function(v1,v2){
    v1minusv2<-v1-v2
    squared.v1minusv2<-v1minusv2*v1minusv2
    out.sqrt<-sqrt(sum(squared.v1minusv2))
    return(out.sqrt)
    
  }### end of fucntion
  
    DNAbin2kmerFreqMatrix<-function(ref,kmer=kmer){ 
    ### return kmer frequency matricies for both ref and que sequences, but only based on kmers found in ref!!!
    ### new kmers in que will be ignored
    #require(seqinr)
    
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
    
    ### do some cleaning by removing IUPAC codes, "-"
    
    mpattern<-"-+[a-z]*"
    
    u.s<-gsub(mpattern,NA,u.s) # strings with "-"
    mpattern<-"[rymkswhbvdn]+"
    u.s<-gsub(mpattern,NA,u.s) # strings with "-"
    
    u.s<- u.s[!is.na(u.s)]
  
    
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
    
    rownames(kmer.freq.matrix)<-rownames(ref)
    if(length(rownames(kmer.freq.matrix))==0) rownames(kmer.freq.matrix)<-names(ref)
    
    out<-list(unique.str=u.s,kmer.Freq.ref=kmer.freq.matrix)
    
    #cat(kmer.freq.matrix)
    return(kmer.freq.matrix)
    
    #return(out)
    
  }#
    
    
  
  #### optimization = FALSE, in this case, a certain kmer will be used, such as kmer = 7
  if (optimization!=TRUE){
    
    #kmer<-7
    kmer.Freq.ref.que<-DNAbin2kmerFreqMatrix2(ref,que,kmer=kmer)
     
    #class(kmer.Freq.ref.que)
    #head(kmer.Freq.ref.que)
    
    #################################################
    ### check model success rate for ref
    #################################################
    set.seed(7)
    knn1<-knn(kmer.Freq.ref.que$kmer.Freq.ref, 
              kmer.Freq.ref.que$kmer.Freq.ref, 
              cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    
    spe.morph<-as.character(Spp2)
    spe.Identified<-as.character(knn1)
    
    spe.morph.Identified<-data.frame(spe.morph,spe.Identified)
    
    matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])
    
    matches<-colSums(matches,dims = 1)
    
    success.rates.ref<-matches[1]/matches[2]
    names(success.rates.ref)<-NULL
    
    #################################################
    ### make prediction for que
    #################################################
    #set.seed(7)
    knn1<-knn(kmer.Freq.ref.que$kmer.Freq.ref, 
              kmer.Freq.ref.que$kmer.Freq.que, 
              cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    
    spe.morph<-as.character(Spp2)
    spe.Identified.que<-as.character(knn1)
    
    #spe.morph.Identified<-data.frame(spe.morph,spe.Identified.que)
    
    #matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])
    
    #matches<-colSums(matches,dims = 1)
    
    #success.rates.que<-matches[1]/matches[2]
    #names(success.rates.que)<-rownames(que)
    #if(length(rownames(que))==0) names(success.rates.que)<-names(que)
    
    #que<-pineMothITS2
    #################################################
    ### calculate theta12 for ref
    
    ref1<-kmer.Freq.ref.que$kmer.Freq.ref    
    rownames(ref1)<-ref.IDs
    if(length(rownames(ref1))==0) names(ref1)<-ref.IDs
    
    FMF.theta12<-FMFtheta12(ref1)
    
    #################################################
    ### calculate FMF for que identified
    
    que1<-kmer.Freq.ref.que$kmer.Freq.que    
    rownames(que1)<-que.IDs
    if(length(rownames(que1))==0) names(que1)<-que.IDs
    
    FMF.que<-numeric(dim(que1)[1])
    
    for(j in 1:dim(que1)[1]){
      #j<-1
      
      #seqInOneSpe<-ref1[grep(spe.Identified[j], FMF.theta12$list.spe, value = FALSE,fixed = TRUE),] 
      
      if(length(rownames(ref1))==0){
        seqInOneSpe<-ref1[grep(spe.Identified.que[j], names(ref1), value = FALSE,fixed = TRUE),]
      }else{seqInOneSpe<-ref1[grep(spe.Identified.que[j], rownames(ref1), value = FALSE,fixed = TRUE),]}
        #seqInOneSpe<-ref1[grep(spe.Identified[j], rownames(ref1), value = FALSE,fixed = TRUE),]
        
      
      
      #class(FMF.theta12$list.spe)
      k<-match(x=spe.Identified.que[j], table=FMF.theta12$list.spe)
      
      #ifelse(FMF.theta12$popSize.PS[k]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) 
      ifelse(FMF.theta12$popSize.PS[k]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) 
      
      
      que2PS.dist<-eucl.dist.two.vect(que1[j,],centroid.spe0)
      #que2PS.dist<-0
      #intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0)
      
      xtheta12<-c(que2PS.dist,FMF.theta12$theta1[k],FMF.theta12$theta2[k])
      FMF.que[j]<-FMF(xtheta12)
    }
    
  
    output.identified<-data.frame(queIDs=que.IDs,
                                  spe.Identified=spe.Identified.que,
                                  FMF=FMF.que)
    out<-list(model.success=success.rates.ref,output_identified=output.identified)
    
    
    #Bayesian.prob=Bayesian.prob)
    #return(output.identified) 
   # return(out) 
    
    
  }else{
    ########################################
    ### kmer.best
    ########################################
    #set.seed(7)
    optimize.kmer<-function (ref,max.kmer=max.kmer){
      #require(ape)
      #set.seed(7)
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
      #seqNames<-NAMES(ref3)
      #ref.IDs<-rownames(ref)
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
        #set.seed(7)
        # i<-1
        
        kmer.freq.ref<-DNAbin2kmerFreqMatrix(ref,kmer=i)
        
        
        knn1<-knn(kmer.freq.ref, kmer.freq.ref, cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
        
        #knn1<-knn1(kmer.freq.ref, kmer.freq.ref, cl=Spp2)
        #attributes(.Last.value)
        #attributes(knn1)
        
        #knn1
        spe.morph<-as.character(Spp2)
        spe.Identified<-as.character(knn1)
        
        spe.morph.Identified<-data.frame(spe.morph,spe.Identified)
        
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
    
    #kmer<-10
    #kmer.optimal<-optimal.kmer(ref,max.kmer=kmer)
    kmer.optimal<-optimize.kmer(ref,max.kmer=kmer)
    success.rates.ref<-kmer.optimal[2]
    kmer.optimal<-kmer.optimal[1]
    
    ########################################
    ### species identification after having got kmer.optimal
    ########################################
    
    kmer.Freq.ref.que<-DNAbin2kmerFreqMatrix2(ref,que,kmer=kmer.optimal)
    
    
    
    #################################################
    ### check model success rate for ref
    #################################################
   # knn1<-knn(kmer.Freq.ref.que$kmer.Freq.ref, 
   #           kmer.Freq.ref.que$kmer.Freq.ref, 
    #          cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    
  
    
   #  spe.morph<-as.character(Spp2)
    # spe.Identified<-as.character(knn1)
    
    # spe.morph.Identified<-data.frame(spe.morph,spe.Identified)
    
    # matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])
    
    # matches<-colSums(matches,dims = 1)
    
    # success.rates.ref<-matches[1]/matches[2]
    # names(success.rates.ref)<-NULL
    
    #################################################
    ### make prediction for que
    #################################################
    #set.seed(7)
    knn1<-knn(kmer.Freq.ref.que$kmer.Freq.ref, 
              kmer.Freq.ref.que$kmer.Freq.que, 
              cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    
    spe.morph<-as.character(Spp2)
    spe.Identified.que<-as.character(knn1)
    
    #spe.morph.Identified<-data.frame(spe.morph,spe.Identified.que)
    
    #matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])
    
    #matches<-colSums(matches,dims = 1)
    
    #success.rates.que<-matches[1]/matches[2]
    #names(success.rates.que)<-rownames(que)
   # if(length(names(success.rates.que))==0) names(success.rates.que)<-names(que)
    #################################################
    ### calculate theta12 for ref
    
    ref1<-kmer.Freq.ref.que$kmer.Freq.ref   
    #class(ref1)
    rownames(ref1)<-ref.IDs
    if(length(rownames(ref1))==0) names(ref1)<-ref.IDs
  
    FMF.theta12<-FMFtheta12(ref1)
    
    #################################################
    ### calculate FMF for que identified
    
    que1<-kmer.Freq.ref.que$kmer.Freq.que    ### matrix
    class(que1)
    rownames(que1)<-que.IDs
    if(length(rownames(que1))==0) names(que1)<-que.IDs
  
    FMF.que<-numeric(dim(que1)[1])
    
    for(j in 1:dim(que1)[1]){
      #j<-1
      
      #seqInOneSpe<-ref1[grep(spe.Identified[j], FMF.theta12$list.spe, value = FALSE,fixed = TRUE),] 
      
      
      if(length(rownames(ref1))==0){
        seqInOneSpe<-ref1[grep(spe.Identified.que[j], names(ref1), value = FALSE,fixed = TRUE),]
      }else{seqInOneSpe<-ref1[grep(spe.Identified.que[j], rownames(ref1), value = FALSE,fixed = TRUE),]}
      
      #seqInOneSpe<-ref1[grep(spe.Identified[j], rownames(ref1), value = FALSE,fixed = TRUE),]
      
      
      #class(FMF.theta12$list.spe)
      k<-match(x=spe.Identified.que[j], table=FMF.theta12$list.spe)
      
      #ifelse(FMF.theta12$popSize.PS[k]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) 
      ifelse(FMF.theta12$popSize.PS[k]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) 
      
      
      que2PS.dist<-eucl.dist.two.vect(que1[j,],centroid.spe0)
      #que2PS.dist<-0
      #intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0)
      
      xtheta12<-c(que2PS.dist,FMF.theta12$theta1[k],FMF.theta12$theta2[k])
      FMF.que[j]<-FMF(xtheta12)
    }
    
    
    
    
    
    output.identified<-data.frame(queIDs=que.IDs,
                                  spe.Identified=spe.Identified.que,
                                  FMF=FMF.que)
    out<-list(model.success=success.rates.ref,output_identified=output.identified)
    #return(output.identified) 
   
    
  }
  
  class(out) <- c("BarcodingR")
  
  return(out) 
}


