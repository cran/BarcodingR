#' Barcoding Gap Calculation
#' 
#' @description Calculation of DNA barcoding gap. Besides K2P distance, raw distance and euclidean could
#' also be used for calculation DNA barcoding gap.
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  dist a character string which takes one of ("raw","K80","euclidean").
#' @return a list indicates the summary statistics of interspecific and intraspecific genetic distance,
#' such as k2P distance.
#' @keywords barcoding gap
#' @export
#' @import ape
#' @import stats
#'
#'
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA, contact at zhangab2008(at)mail.cnu.edu.cn
#' @note the current version of the function can only be used for protein-coding barcodes,
#' such as, COI. The futuren version may incorporate calculation for non-coding barcodes,for
#' instance, ITS1, ITS2.
#' @references
#' Meyer, Christopher P., and Gustav Paulay. (2005). ''DNA barcoding: error rates based on
#' comprehensive sampling.''.PLoS biology 3.12: e422.
#'
#' F.Jiang, Q. Jin, L. Liang, A.B. Zhang,and Z.H. Li.(2014). Existence of Species
#' Complex Largely Reduced Barcoding Success for Invasive Species of Tephritidae: A Case Study in
#' Bactrocera spp. Mol Ecol Resour. 14(6):1114-1128 DOI: 10.1111/1755-0998.12259.
#'
#'
#' @examples
#'
#' data(TibetanMoth)
#' TibetanMoth<-as.DNAbin(as.character(TibetanMoth[1:20,]))
#' b.gap<-barcoding.gap(ref=TibetanMoth,dist="K80")
#' b.gap


barcoding.gap<-function (ref,dist=dist){ ### raw/k2p,"raw" "k2p", "euclidean"

  maxInDist <-function(distobj, sppVector = NULL, propZero = FALSE, rmNA = FALSE){
      dat <- as.matrix(distobj)
      if(length(sppVector) > 0) dimnames(dat)[[1]] <- sppVector
      conSpecDists <- list()
      for (i in 1:length(dimnames(dat)[[1]])) {
        conSpec <- dimnames(dat)[[1]] == dimnames(dat)[[1]][i]
        conSpecDists[[i]] <- max(dat[conSpec, i], na.rm = rmNA)
      }
      if (propZero)
        output <- length(which(unlist(conSpecDists) == 0))/length(unlist(conSpecDists))
      else output <- unlist(conSpecDists)
      output
    }
  nonConDist <-function(distobj, sppVector = NULL, propZero = FALSE, rmNA = FALSE){
      distobj <- as.matrix(distobj)
      if(length(sppVector) > 0) dimnames(distobj)[[1]] <- sppVector
      nonSpecDists <- list()
      for(i in 1:length(dimnames(distobj)[[1]])){
        nonSpec <- dimnames(distobj)[[1]] != dimnames(distobj)[[1]][i]
        nonSpecDists[[i]] <- min(distobj[nonSpec,i] , na.rm = rmNA)
      }
      if(propZero) output <- length(which(unlist(nonSpecDists) == 0))/length(unlist(nonSpecDists)) else output <- unlist(nonSpecDists)

      output
    }


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
   twoSpeDist<-function(sp1,sp2,dist=dist){ ### "raw" "k2p", "euclidean"
    ### sp1,sp2, DNAbin
    ###

    n1<-dim(sp1)[1]
    n2<-dim(sp2)[1]

    sp12<-rbind(sp1,sp2)

    if(dist=="euclidean"){

      digitize.DNA<-function(seqs){

        locus<-toupper(as.character(seqs))
        digitized.DNA<-locus
        digitized.DNA[digitized.DNA=="A"]<-0.1
        digitized.DNA[digitized.DNA=="T"]<-0.2
        digitized.DNA[digitized.DNA=="G"]<-0.3
        digitized.DNA[digitized.DNA=="C"]<-0.4
        digitized.DNA[digitized.DNA=="-"]<-0.5
        digitized.DNA[digitized.DNA=="N"]<-0
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


      sp12.digitized<-digitize.DNA(sp12)

      #dist<-dist(sp12.digitized, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

      dist<-dist(sp12.digitized, method = "euclidean", diag = F, upper = F, p = 2)

      #dist<-dist.dna(sp12,model = "raw")  ###  "raw" "K80"
      #dist
      # class(dist)

      dist<-as.matrix(dist)
      diag(dist) <- NA

      #lower.tri(dist,diag = FALSE)

      inter<-dist[(n1+1):(n1+n2),1:n1]
      #inter
      inter<-list(inter)
      inter<-unlist(inter)
      #inter


      intra1<-dist[1:n1,1:n1]
      intra1[upper.tri(intra1)] <- NA
      #intra1

      intra1<-list(intra1)
      intra1<-unlist(intra1)
      #intra1

      intra2<-dist[(n1+1):(n1+n2),(n1+1):(n1+n2)]
      intra2[upper.tri(intra2)] <- NA
      # intra2

      intra2<-list(intra2)
      intra2<-unlist(intra2)

      intra12<-c(intra1,intra2)
      intra12<-intra12[!is.na(intra12)]

      #intra12
      out<-list(intra12,inter)


    }else{ ### "raw"  "k2p"

      #dist<-dist.dna(sp12,model = "K80")  ###  "raw" "K80"
      dist<-dist.dna(sp12,model = dist)  ###  "raw" "K80"
      dist<-as.matrix(dist)
      diag(dist) <- NA
      inter<-dist[(n1+1):(n1+n2),1:n1]
      inter<-list(inter)
      inter<-unlist(inter)
      intra1<-dist[1:n1,1:n1]
      intra1[upper.tri(intra1)] <- NA
      #intra1

      intra1<-list(intra1)
      intra1<-unlist(intra1)
      #intra1

      intra2<-dist[(n1+1):(n1+n2),(n1+1):(n1+n2)]
      intra2[upper.tri(intra2)] <- NA
      # intra2

      intra2<-list(intra2)
      intra2<-unlist(intra2)

      intra12<-c(intra1,intra2)
      intra12<-intra12[!is.na(intra12)]

      #intra12
      out<-list(intra12,inter)
    }
    return(out)
  } ### the end of function twoSpeDist() ### "euclidean"
  multSpeDist<-function(ref,dist=dist){     ###iple
    sampleSpeNames<-attr(ref,"dimnames")[[1]]
    mpattern<-".+,"
    #mpattern<-".+,Noctuidae_"
    #mpattern<-"Noctuidae_"
    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    #Spp
    f<-factor(Spp)
    intra<-NULL
    inter<-NULL

    for(i in 1:(length(levels(f))-1)){
      cat(paste("i=",i),"\n")
      seqInOneSpe<-ref[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),]

      for(j in (i+1):length(levels(f))){
        seqInOneSpe2<-ref[grep(levels(f)[j], Spp, value = FALSE,fixed = TRUE),]

        #t.d<-twoSpeDist(seqInOneSpe,seqInOneSpe2,dist="raw") ### "euclidean","K80",""
        t.d<-twoSpeDist(seqInOneSpe,seqInOneSpe2,dist=dist) ### "euclidean","K80",""

        intra<-c(intra,t.d[[1]])
        inter<-c(inter,t.d[[2]])
      }
    }

    out2<-list(intra,inter)
    return(out2)

  } ### "raw"

  m.d<-multSpeDist(ref,dist=dist)  ###  "raw", "K80"

  ref.IDs<-NAMES(ref)

  morph.spe<-gsub(".+,","",ref.IDs) # remove sequence ID before ","
  k2p.dist<-dist.dna(ref,model = "raw",pairwise.deletion = TRUE)
  inter<-nonConDist(k2p.dist,morph.spe,rmNA = T)
  intra<-maxInDist(k2p.dist,morph.spe,rmNA = T)
  intra<-m.d[[1]]
  inter<-m.d[[2]]
  h.inter<-hist(inter,freq = TRUE)
  h.intra<-hist(intra,freq = TRUE)
  ### plot
  max.x<-max(c(inter,intra))
  #max.x<-median(c(inter,intra))*2
  max.y<-max(c(h.inter$counts,h.intra$counts))
   # plot(c(0,max.x),c(0,max.y*0.1),type="n", ## ver1.02
     plot(c(0,max.x),c(0,max.y*0.5),type="n",
       xlab="genetic distance",ylab="Frequency",
       main="DNA barcoding gap analysis",
       #sub="red-intra,blue-inter"
       )
  #title("DNA barcoding gap analysis")

  h.inter2<-hist(inter,freq = TRUE,breaks = "Sturges",col = "blue",border = "white",add = TRUE)
  #h.inter2<-hist(inter,freq = TRUE,breaks = 100,col = "blue",border = "white",add = TRUE)## ver1.02
  #h.inter2<-hist(inter,freq = TRUE,col = "blue",border = "white",add = TRUE)
  #h.inter2<-hist(inter,freq = TRUE,breaks = 12,col = "blue",border = "white")
  h.inter.xfit<-seq(min(inter),max(inter),length=40)
  h.inter.yfit<-dnorm(h.inter.xfit,mean=mean(h.inter.xfit),sd=sd(h.inter.xfit))
  h.inter.yfit<-h.inter.yfit*diff(h.inter2$mids[1:2])*length(inter)
  lines(h.inter.xfit,h.inter.yfit,col="blue",lwd=2)


  h.intra2<- hist(intra,freq = TRUE,breaks = "Sturges",col = "red",border = "white",add = TRUE)
  #h.intra2<- hist(intra,freq = TRUE,breaks = 100,col = "red",border = "white",add = TRUE) ## ver1.02
  #h.intra2<- hist(intra,freq = TRUE,col = "red",border = "white",add = TRUE)

  h.intra.xfit<-seq(min(intra),max(intra),length=40)
  h.intra.yfit<-dnorm(h.intra.xfit,mean=mean(h.intra.xfit),sd=sd(h.intra.xfit))
  h.intra.yfit<-h.intra.yfit*diff(h.intra2$mids[1:2])*length(intra)
  lines(h.intra.xfit,h.intra.yfit,col="red",lwd=2)

  legend(x = 0.65*max.x,y = 0.5*max.y,
         #legend(x = 0.70*max.x,y = 0.10*max.y, ## ver1.02
         legend = c("intraspecific","interspecific"),
         pch = c(22,22),
         col = c("red","blue"),
         #bg = c("red","blue")
         )

  ##################

  out<-list(c("inter","intra"),c(summary(inter),summary(intra)))

  return(out)

}


#b.gap<-barcoding.gap(ref,dist="raw")
#b.gap<-barcoding.gap(ref,dist="K80")
#b.gap<-barcoding.gap(ref,dist="euclidean")































