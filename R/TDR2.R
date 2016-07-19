#' TDR2 Species Membership Value 
#' 
#' @description To calculate TDR value for a set of queries and one potential species. Its value is in the range of [0,1],
#' 0 indicates extremly weak species membership, values close 1 indicating strong species membership.
#' 
#' @param  oneSpe object of class "DNAbin" which contains DNA squences from one speices
#' @param  que object of class "DNAbin" which contains DNA squences different samples
#' @param  boot a numeric value indicating times of resampling along sequence columns
#' @param  boot2 a numeric value indicating times of resampling along sequence rows (different samples)
#' @return a numeric vector represents TDR values for each query against the species
#' @keywords TDR2
#' @export 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @references Jin Q, L,J.He, A.B. Zhang* (2012). A Simple 2D Non-Parametric Resampling Statistical Approach to 
#' Assess Confidence in Species Identification in DNA Barcoding-An Alternative to Likelihood and Bayesian Approaches.
#' PLoS ONE 7(12): e50831. doi:10.1371/journal.pone.0050831.http://dx.plos.org/10.1371/journal.pone.0050831.
#' @note oneSpe and que should be the same in sequence length, i.e., they should be aligned in prior.
#' It's strongly recommended that oneSpe should have large enough sample size,e.g., 20.
#'  
#' @examples
#' 
#' data(TibetanMoth)  
#' sampleSpeNames<-NAMES(TibetanMoth)
#' Spp<-gsub(".+,","",sampleSpeNames)
#' oneSpe<-TibetanMoth[grep("Macdunnoughia_crassisigna", Spp, value = FALSE,fixed = TRUE),] 
#' oneSpe<-as.DNAbin(as.character(oneSpe[1:5,]))
#' que<-TibetanMoth[grep("Agrotis_justa", Spp, value = FALSE,fixed = TRUE),] 
#' que2<-oneSpe[1:2,]
#' out<-TDR2(oneSpe,que, boot=10,boot2=10) ### true false identification



TDR2<-function (oneSpe,que, boot,boot2){
  
  
  if (class(oneSpe)!="DNAbin"||class(que)!="DNAbin") 
    stop("invalid sequence format! DNAbin format is required for DNA seqs!")
  
  if (dim(que)[2] != dim(oneSpe)[2]) 
    warning("sequences in ref and que are different in length!")
  
  no.indi<-dim(oneSpe)[1]
  no.que<-dim(que)[1]
  
  if (no.indi<1) 
    stop("oneSpe has at least 2 individuals, 20 recommended!")
  
  
  
  
  
  digitize.DNA<-function(seqs){
    
    locus<-toupper(as.character(seqs))
    digitized.DNA<-locus
    digitized.DNA[digitized.DNA=="A"]<-0.1
    digitized.DNA[digitized.DNA=="T"]<-0.2
    digitized.DNA[digitized.DNA=="G"]<-0.3
    digitized.DNA[digitized.DNA=="C"]<-0.4
    digitized.DNA[digitized.DNA=="-"]<-0.5
    digitized.DNA[digitized.DNA=="N"]<-0.6
    digitized.DNA[digitized.DNA=="R"]<-0
    digitized.DNA[digitized.DNA=="Y"]<-0
    digitized.DNA[digitized.DNA=="M"]<-0
    digitized.DNA[digitized.DNA=="K"]<-0
    digitized.DNA[digitized.DNA=="S"]<-0
    digitized.DNA[digitized.DNA=="W"]<-0
    digitized.DNA[digitized.DNA=="H"]<-0
    digitized.DNA[digitized.DNA=="B"]<-0
    digitized.DNA[digitized.DNA=="V"]<-0
    digitized.DNA[digitized.DNA=="D"]<-0
    
    
    digitized.DNA<-as.numeric(digitized.DNA)
    #digitized.DNA<-as.matrix(digitized.DNA)
    digitized.DNA2<-array(digitized.DNA,dim=dim(seqs))
    dim(digitized.DNA2)
    
    return(digitized.DNA2)
  }
  eucl.dist.two.vect<-function(v1,v2){
    v1minusv2<-v1-v2
    squared.v1minusv2<-v1minusv2*v1minusv2
    out.sqrt<-sqrt(sum(squared.v1minusv2))
    return(out.sqrt)
    
  }### end of fucntion
  
  
  
  
  
  
  locus.as.character<-as.character(oneSpe)
  locus.as.character<-array(locus.as.character,c(dim(oneSpe)[1],dim(oneSpe)[2]))
  
  locus.as.character2<-as.character(que)
  locus.as.character2<-array(locus.as.character2,c(dim(que)[1],dim(que)[2]))
  
  
  
  
  #ks.best<-numeric(boot)
  #sil.ks.best<-numeric(boot)
  
  ### 
  
  TDR<-numeric(no.que)
  
  for (j in 1:no.que){
    # j<-1
    oneSpe.que<-rbind(locus.as.character,locus.as.character2[j,]) ### oneSpe + one que
    
    accepted.times<-0
    for (i in 1:boot) { ### i for oneSpe
    # i<-1
      cat("i:",i,"\n")
      
    ### generate each resampled DNA sequence matrix for the species
    #resampled.locus.as.character<-locus.as.character[,sample(ncol(locus.as.character), replace = TRUE)]
    resampled.locus.as.character<-oneSpe.que[,sample(dim(oneSpe)[2], replace = TRUE)]
    
    #
      # j<-1
     # resampled.locus.as.character.plus.one.query<-rbind(resampled.locus.as.character,
       #                                                  locus.as.character2[j,]
        #                                                 )
      #dim(resampled.locus.as.character.plus.one.query)
      #resampled.locus.as.character.plus.one.query2<-as.DNAbin(resampled.locus.as.character.plus.one.query)
      resampled.locus.as.character.plus.one.query2<-as.DNAbin(resampled.locus.as.character)
    
       
    
      #rlacpoq.digitized
      
      ### calculate the distance of the que to the species firstly!
      rlacpoq.digitized<-digitize.DNA(resampled.locus.as.character.plus.one.query2)
      
      seqInOneSpe<-rlacpoq.digitized[1:no.indi,]
      que.tmp<-rlacpoq.digitized[(no.indi+1),]
      
      #centroid.spe<-apply(seqInOneSpe,2,mean)
      
      inter.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=que.tmp)
      
      dist.que2spe<-min(inter.dist)
      
      
      #dist.que2spe<-eucl.dist.two.vect(que.tmp,centroid.spe)
      
      cat("dist.que2spe:",dist.que2spe,"\n")
      
      
      ### vertical sampling for boot2 times with replacement
        
      dist.simu<-numeric(boot2)  
    
      for (k in 1:boot2){
        
      resampled2.locus.as.character.plus.one.query<-resampled.locus.as.character[sample(dim(oneSpe)[1]+1, replace = TRUE),]
      #resampled2.locus.as.character.plus.one.query<-resampled.locus.as.character.plus.one.query[sample(dim(oneSpe)[1]+1, replace = F),]
      
      
      
      
      resampled2.locus.as.character.plus.one.query2<-as.DNAbin(resampled2.locus.as.character.plus.one.query)
      #rlacpoq.digitized
      
      ### calculate the distance of the que to the species firstly!
      rlacpoq.digitized<-digitize.DNA(resampled2.locus.as.character.plus.one.query2)
      
      seqInOneSpe<-rlacpoq.digitized[1:no.indi,]
      que.tmp2<-rlacpoq.digitized[(no.indi+1),]
      
      #centroid.spe<-apply(seqInOneSpe,2,mean)
      
      inter.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=que.tmp2)
      
      dist.simu[k]<-min(inter.dist)
      
      
      #dist.simu[k]<-eucl.dist.two.vect(que.tmp2,centroid.spe)
      
      
      } ### the end of k-loop
      
      #ifelse(popSize.PS[1]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean)) ###
      ci_dist.simu<-quantile(dist.simu,prob=c(0.025,0.975))
      #ci_dist.simu[1]  
      
      cat("ci_dist.simu",ci_dist.simu,"\n")
    
      ifelse(ci_dist.simu[1]<=dist.que2spe&&dist.que2spe<=ci_dist.simu[2],
             accepted.times<-accepted.times+1,
             accepted.times<-accepted.times)
      
      dist.simu
      
      
    }### the end of i-loop
    
    
    
    
    TDR[j]<-accepted.times/boot
    cat(j,":\n")
    cat("accepted.times",accepted.times,"\n")
    
    
  } ### the end of j-loop
  
  return(TDR)
  
}











