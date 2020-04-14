#' Calculate Intraspecific and Interspecific Variation
#' 
#' @description Calculation intraspecific variation (sd) of the potential species theta1, and mean interspecific
#' distance (here, the mean distance between the potential species and its nearest neighbor theta2) 
#' (fuzzy-set based method,slightly modified from Zhang et al. 2012). The calculation was done for all species in the 
#' reference dataset.
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @return a data frame containing intraspecific (sd, theta1) and interspefic variation (mean) of all species,
#' and their corresponding nearest neighbor (NN).
#' @keywords FMFtheta12
#' @export 
#' @import stats
#' @import graphics
#' @import utils 
#' @import sp
#' 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA, contact at zhangab2008 (at) mail.cnu.edu.cn.
#' @references 
#' Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach 
#' to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
#' @examples
#' 
#' data(TibetanMoth) 
#' ref<-as.DNAbin(as.character(TibetanMoth[1:50,]))
#' FMF.theta12<-FMFtheta12(ref)
#' FMF.theta12

FMFtheta12<-function(ref){
  
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
  
  ### 2.1 dealing with input DNA data! 
  ### 2.2 calculate species center vectors  
  #morph.spe<-gsub(".+,","",rownames(Ref)) # remove sequence ID before "," 
  morph.spe<-gsub(".+,","",NAMES(ref)) # remove sequence ID before ","
  #no.morph.spe<-length(unique(morph.spe))
  
  ref<-ref[,seg.sites(ref)]
  ref<-digitize.DNA(ref)
  rownames(ref)<-morph.spe #sampleSpeNames
  
  
  #species.centers<-aggregate(scale(Ref),by=list(morph.spe),FUN="mean")
  species.centers<-aggregate(ref,by=list(morph.spe),FUN="mean")
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
  uniSpeNames<-unique(Spp)
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
  
  seqInOneSpe<-ref[grep(levels(f)[1], Spp, value = FALSE,fixed = TRUE),]
  #dim(seqInOneSpe)==NULL
  ifelse(popSize.PS[1]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean)) ###
  ifelse(popSize.PS[1]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) ###  centroid.spe0 - 
  
  length.sites<-length(centroid.spe)
  ifelse(popSize.PS[1]==1,intra.dist<-0,intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe)) ### to the centroid
  ifelse(popSize.PS[1]==1,theta1.tmp[1]<-0,theta1.tmp[1]<-3*sd(intra.dist)) ### theta1 is sligthly different from 
  #ifelse(popSize.PS[1]==1,theta1.tmp[1]<-0,theta1.tmp[1]<-max(intra.dist)) ### theta1 is sligthly different from 
    for(i in 2:length(levels(f))){
    cat(paste("i=",i),"\n")
   
    #seqInOneSpe<-sDNAbin[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),]
    seqInOneSpe<-ref[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),] 
    ifelse(popSize.PS[i]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) 
    #ifelse(popSize.PS[i]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean))
    ifelse(popSize.PS[i]==1,centroid.spe<-c(centroid.spe,seqInOneSpe),centroid.spe<-c(centroid.spe,apply(seqInOneSpe,2,mean)))
    #ifelse(dim(seqInOneSpe)==NULL,theta1.tmp[i]<-0,theta1.tmp[i]<-max(dist(seqInOneSpe))) # error!
    #ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-max(dist(seqInOneSpe)))
    
    ifelse(popSize.PS[i]==1,intra.dist<-eucl.dist.two.vect(seqInOneSpe,centroid.spe0),intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0))
    
    
    #intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0) ### to the centroid 
    #ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-max(intra.dist)) ### theta1 is sligthly different from 
    ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-3*sd(intra.dist)) ### theta1 is sligthly different from 
    
  }
  centroid.spe.matrix<-t(array(centroid.spe,dim=c(length.sites,length(centroid.spe)%/%length.sites)))
  
  ### 1.2 calculate theta2 for each pairt of species
  ###################
  #codes<-out.somu$out.som.unique$codes  ### seqsRef$codes
  codes<-centroid.spe.matrix
   theta12.all<-data.frame(list.spe=list.spe,PS=pairs[,1],NN=pairs[,2],stringsAsFactors=TRUE)
  #source("eucl.dist.two.vect.R")
  theta2.tmp<-numeric(dim(pairs)[1])
  
  
  for (i in 1:dim(pairs)[1]){
    v1<-codes[pairs[i,1],]
    v2<-codes[pairs[i,2],]
    theta2.tmp[i]<-eucl.dist.two.vect(v1,v2)
  }
  
  theta12.all$theta1<-with(theta12.all,theta1.tmp)
  theta12.all$theta2<-with(theta12.all,theta2.tmp)
  
  theta12.all$popSize.PS<-with(theta12.all,popSize.PS)
  return(theta12.all) 
  
} ### the end of the function
























