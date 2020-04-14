

#' Species Identification using Protein-coding Barcodes
#'
#' @description Species identification using protein-coding barcodes with different methods,including BP-based method
#' (Zhang et al. 2008), fuzzy-set based method (Zhang et al. 2012), Bayesian-based method (Jin et al. 2013).
#'
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  que object of class "DNAbin", whose identities (species names) need to be inferred.
#' @param  method a character string indicating which method will be used to train model and/or infer species membership. 
#'         One of these methods ("fuzzyId", "bpNewTraining", "bpNewTrainingOnly", "bpUseTrained","Bayesian") should be specified.
#' 
#' @return a list containing model parameters used, species identification success rates using references,
#'         query sequences, species inferred, and corresponding confidence levels 
#'         (bp probability for BP-based method / FMF values for fuzzy set theory based method / posterior probability for Bayesian method) when available.
#' 
#' @keywords BSI
#' 
#' @export
#' 
#' @import ape
#' @import nnet
#' @import class
#' @import stats
#' @import sp
#'
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA.
#'         zhangab2008(at)mail.cnu.edu.cn
#' 
#' @note functions fasta2DNAbin() from package:adegenet and read.dna() from package:ape were used to obtain DNAbin object in our package.
#' The former is used to read large aligned coding DNA barcodes, the latter unaligned ones. ref and que
#' should be aligned with identical sequence length. We provided a pipeline to perform fast
#' sequences alignment for reference and query sequences. Windows users could contact zhangab2008(at)mail.cnu.edu.cn
#' for an exec version of the package. For very large DNA dataset, read.fas() package:phyloch is strongly suggested instead of
#' fasta2DNAbin() since the latter is very slow.
#' 
#' @references
#' Zhang, A. B., M. D. Hao, C. Q. Yang, and Z. Y. Shi. (2017). BarcodingR: an integrated R package for species identification using DNA barcodes. Methods Ecol Evol. 8(5):627-634.
#' https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12682.
#' 
#' Jin,Q., H.L. Han, X.M. Hu, X.H. Li,C.D. Zhu,S. Y. W. Ho, R. D. Ward, A.B. Zhang . (2013). Quantifying Species Diversity with a DNA Barcoding-Based Method: Tibetan Moth Species (Noctuidae) on the Qinghai-Tibetan Plateau. PloS One 8: e644.
#' https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064428.
#' 
#' Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
#' https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2011.05235.x
#' 
#' Zhang, A. B., D. S. Sikes, C. Muster, S. Q. Li. (2008). Inferring Species Membership using DNA sequences with Back-propagation Neural Networks. Systematic Biology, 57(2):202-215.
#' https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12682

#' @examples
#' data(TibetanMoth)
#' ref<-as.DNAbin(as.character(TibetanMoth[1:5,]))
#' que<-as.DNAbin(as.character(TibetanMoth[50:55,]))
#' bsi<-barcoding.spe.identify(ref, que, method = "fuzzyId")
#' bsi
#' bsi<-barcoding.spe.identify(ref, que, method = "bpNewTraining")
#' bsi
#' bsi<-barcoding.spe.identify(ref, que, method = "Bayesian")
#' bsi



#library(BarcodingR)

barcoding.spe.identify<-function(ref, que, method = "bpNewTraining") {##"bpNewTraining", "bpUseTrained", "fuzzyId","Bayesian"
  #barcoding.species.identification<-function (ref, que, method = "bp") {##"fuzzyId","bp","Bayesian"
 # "bpNewTraining", ### train and identify at the same time
 # "bpNewTrainingOnly" ### just train the model for later use
 # "bpUseTrained", ### use the trained model in "bpNewTrainingOnly" which is save to a tmp file.

  if (dim(que)[2] != dim(ref)[2])
    warning("sequences in ref and que are different in length!")

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

  strings.equal<-function(str1,str2){ifelse(str1==str2,1,0)}

  ##############
  naiveBayes <- function(x, ...)
    UseMethod("naiveBayes")

  naiveBayes.default <- function(x, y, laplace = 0, ...) {
    call <- match.call()
    Yname <- deparse(substitute(y))
    x <- as.data.frame(x,stringsAsFactors=TRUE)

    ## estimation-function
    est <- function(var)
      if (is.numeric(var)) {
        cbind(tapply(var, y, mean, na.rm = TRUE),
              tapply(var, y, sd, na.rm = TRUE))
      } else {
        tab <- table(y, var)
        (tab + laplace) / (rowSums(tab) + laplace * nlevels(var))
      }

    ## create tables
    apriori <- table(y)
    tables <- lapply(x, est)

    ## fix dimname names
    for (i in 1:length(tables))
      names(dimnames(tables[[i]])) <- c(Yname, colnames(x)[i])
    names(dimnames(apriori)) <- Yname

    structure(list(apriori = apriori,
                   tables = tables,
                   levels = levels(y),
                   call   = call
    ),

    class = "naiveBayes"
    )
  }

  naiveBayes.formula <- function(formula, data, laplace = 0, ...,
                                 subset, na.action = na.pass) {
    call <- match.call()
    Yname <- as.character(formula[[2]])

    if (is.data.frame(data)) {
      ## handle formula
      m <- match.call(expand.dots = FALSE)
      m$... <- NULL
      m$laplace = NULL
      m$na.action <- na.action
      m[[1]] <- as.name("model.frame")
      m <- eval(m, parent.frame())
      Terms <- attr(m, "terms")
      if (any(attr(Terms, "order") > 1))
        stop("naiveBayes cannot handle interaction terms")
      Y <- model.extract(m, "response")
      X <- m[,-attr(Terms, "response"), drop = FALSE]

      return(naiveBayes(X, Y, laplace = laplace, ...))
    } else if (is.array(data)) {
      nam <- names(dimnames(data))
      ## Find Class dimension
      Yind <- which(nam == Yname)

      ## Create Variable index
      deps <- strsplit(as.character(formula)[3], ".[+].")[[1]]
      if (length(deps) == 1 && deps == ".")
        deps <- nam[-Yind]
      Vind <- which(nam %in% deps)

      ## create tables
      apriori <- margin.table(data, Yind)
      tables <- lapply(Vind,
                       function(i) (margin.table(data, c(Yind, i)) + laplace) /
                         (as.numeric(apriori) + laplace * dim(data)[i]))
      names(tables) <- nam[Vind]

      structure(list(apriori = apriori,
                     tables = tables,
                     levels = names(apriori),
                     call   = call
      ),

      class = "naiveBayes"
      )
    } else stop("naiveBayes formula interface handles data frames or arrays only")

  }


  print.naiveBayes <- function(x, ...) {
    cat("\nNaive Bayes Classifier for Discrete Predictors\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\nA-priori probabilities:\n")
    print(x$apriori / sum(x$apriori))

    cat("\nConditional probabilities:\n")
    for (i in x$tables) {print(i); cat("\n")}

  }

  predict.naiveBayes <- function(object,
                                 newdata,
                                 type = c("class", "raw"),
                                 threshold = 0.001,
                                 eps = 0,
                                 ...) {
    type <- match.arg(type)
    newdata <- as.data.frame(newdata,stringsAsFactors=TRUE)
    attribs <- match(names(object$tables), names(newdata))
    isnumeric <- sapply(newdata, is.numeric)
    newdata <- data.matrix(newdata)
    L <- sapply(1:nrow(newdata), function(i) {
      ndata <- newdata[i, ]
      L <- log(object$apriori) + apply(log(sapply(seq_along(attribs),
                                                  function(v) {
                                                    nd <- ndata[attribs[v]]
                                                    if (is.na(nd)) rep(1, length(object$apriori)) else {
                                                      prob <- if (isnumeric[attribs[v]]) {
                                                        msd <- object$tables[[v]]
                                                        msd[, 2][msd[, 2] <= eps] <- threshold
                                                        dnorm(nd, msd[, 1], msd[, 2])
                                                      } else object$tables[[v]][, nd]
                                                      prob[prob <= eps] <- threshold
                                                      prob
                                                    }
                                                  })), 1, sum)
      if (type == "class")
        L
      else {
        ## Numerically unstable:
        ##            L <- exp(L)
        ##            L / sum(L)
        ## instead, we use:
        sapply(L, function(lp) {
          1/sum(exp(L - lp))
        })
      }
    })
    if (type == "class")
      factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
    else t(L)
  }



  ##############





  #source("aggregate.R")


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

  #[1] 0.6289163 0.1465522 0.6379375


  sampleSpeNames<-attr(ref,"dimnames")[[1]]

  mpattern<-".+,"
  #mpattern<-".+,Noctuidae_"
  #mpattern<-"Noctuidae_"
  #mpattern<-".+,[[:alpha:]]+_"

  Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
  #Spp



  bp<-function(ref1,sampleSpeNames2, method = "bpNewTraining",que1){
   # "bpNewTraining", ### train and identify at the same time
   # "bpNewTrainingOnly" ### just train the model for later use
   # "bpUseTrained",

    if(method == "bpNewTraining"){

    ref1_tmp<-ref1[!is.na(ref1)]
    range1<-1./max(abs(ref1_tmp))
    n.hidden<-ceiling(log2(dim(ref1)[1]))
    nnet.trained <- nnet(ref1, sampleSpeNames2, size = n.hidden, rang = range1,
              decay = 5e-5,
              maxit = 1e+6,
              abstol = 1.0e-8,
              reltol = 1.0e-8,
              MaxNWts = 20000)
              #decay = 5e-4, maxit = 1e+5,MaxNWts = 2000)

  spe.inferred0<-predict(nnet.trained, que1)
  spe.inferred<-spe.inferred0
  #spe.inferred[spe.inferred>=0.95]<-1
  #spe.inferred[spe.inferred<0.95]<-0
  #spe.inferred
  #colnames(spe.inferred)
  #which.max(spe.inferred[1,])
  inferred<-apply(spe.inferred,1,FUN=which.max)
  inferred.prob<-apply(spe.inferred,1,FUN=max)
  #inferred
  colnames(spe.inferred)[inferred]
  #spe.inferred0<-t(spe.inferred0)
  #bp.prob<-spe.inferred0[inferred]
  output.identified<-data.frame(queIDs=queIDs,
                              spe.Identified=colnames(spe.inferred)[inferred],
                              bp.prob=inferred.prob,stringsAsFactors=TRUE)


    rownames(output.identified)<-NULL

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

    spe.morph.Identified<-data.frame(Spp,p.ref,stringsAsFactors=TRUE)

    matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])


    matches<-colSums(matches,dims = 1)

    success.rates.ref<-matches[1]/matches[2]
    names(success.rates.ref)<-NULL

    ############################################################
    ####### calculate model success rate using ref the end.
    ###########################################################

    out.bp<-list(summary.model=summary(nnet.trained),###convergence=nnet.trained$convergence,
                 convergence=nnet.trained$convergence,
                 success.rates.ref=success.rates.ref,
                 output_identified=output.identified)



    #current.wd<-getwd()
    #setwd
    
    #Rhome<-R.home() ### 2020/4/13 21:06:46
    Rhome<-tempdir() ### 2020/4/14 21:31:02
    fileName<-"bbsi_tmp"
    fileName<-paste(Rhome,fileName,sep = "/")### 2020/4/13 21:06:46
    fileName
    #fileName<-paste("simulation",i,sep = "")
    fileName<-paste(fileName,".RData",sep = "")
    fileName
    
    save(nnet.trained,
         out.bp,
         success.rates.ref,
         #unique.str.ref,
         #center.ref1,
         file = fileName) ## #save(x, y, file = "xy.RData")

    #return(out.bp)
    output2.identified<-out.bp

    }
    if(method == "bpNewTrainingOnly"){



      ref1_tmp<-ref1[!is.na(ref1)]
      center.ref1<-apply(ref1,MARGIN=2,FUN=mean)

      range1<-1./max(abs(ref1_tmp))
      n.hidden<-ceiling(log2(dim(ref1)[1]))
      nnet.trained <- nnet(ref1, sampleSpeNames2, size = n.hidden, rang = range1,
                           decay = 5e-5,
                           maxit = 1e+6,
                           abstol = 1.0e-8,
                           reltol = 1.0e-8,
                           MaxNWts = 20000)

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

      spe.morph.Identified<-data.frame(Spp,p.ref,stringsAsFactors=TRUE)

      matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])


      matches<-colSums(matches,dims = 1)

      success.rates.ref<-matches[1]/matches[2]
      names(success.rates.ref)<-NULL

      ############################################################
      ####### calculate model success rate using ref the end.
      ###########################################################

     # fileName<-"bbsi_tmp"

      #fileName<-paste(fileName,".RData",sep = "")
      
     #Rhome<-R.home() ### 2020/4/13 21:06:46
      Rhome<-tempdir() ### 2020/4/14 21:31:02
      fileName<-"bbsi_tmp"
      fileName<-paste(Rhome,fileName,sep = "/")### 2020/4/13 21:06:46
      fileName
      #fileName<-paste("simulation",i,sep = "")
      fileName<-paste(fileName,".RData",sep = "")
      fileName
      
      
      
      
      save(nnet.trained,
           success.rates.ref,
           #center.ref1,
           file = fileName) ## #save(x, y, file = "xy.RData")

      output2.identified<-"just saved to file!"


    }###
    if(method == "bpUseTrained"){

      ##################### ### 2020/4/13 21:06:46
      #Rhome<-R.home() ### 2020/4/13 21:06:46
      Rhome<-tempdir() ### 2020/4/14 21:31:02
      fileName<-"bbsi_tmp"
      fileName<-paste(Rhome,fileName,sep = "/")### 2020/4/13 21:06:46
      fileName
      #fileName<-paste("simulation",i,sep = "")
      fileName<-paste(fileName,".RData",sep = "")
      fileName
      ############################ 2020/4/13 21:06:46
      
      #if(!file.exists("bbsi_tmp.RData"))
        if(!file.exists(fileName)) ###
        stop("file bbsi_tmp.RData not found in R.home directory!
         you need to rebuild the model!
         ")


      #load("bbsi_tmp.RData")
      load(fileName) ####2020/4/13 21:14:15
      
      queIDs<-attr(que,"dimnames")[[1]]



      spe.inferred0<-predict(nnet.trained, que1)
      spe.inferred<-spe.inferred0
      #spe.inferred[spe.inferred>=0.95]<-1
      #spe.inferred[spe.inferred<0.95]<-0
      #spe.inferred
      #colnames(spe.inferred)
      #which.max(spe.inferred[1,])
      inferred<-apply(spe.inferred,1,FUN=which.max)
      inferred.prob<-apply(spe.inferred,1,FUN=max)
      #inferred
      colnames(spe.inferred)[inferred]
      #spe.inferred0<-t(spe.inferred0)
      #bp.prob<-spe.inferred0[inferred]
      output.identified<-data.frame(queIDs=queIDs,
                                    spe.Identified=colnames(spe.inferred)[inferred],
                                    bp.prob=inferred.prob,stringsAsFactors=TRUE)


      rownames(output.identified)<-NULL

      output2.identified<-list(summary.model=summary(nnet.trained),###convergence=nnet.trained$convergence,
                   convergence=nnet.trained$convergence,
                   success.rates.ref=success.rates.ref,
                   output_identified=output.identified)



    }



    class(output2.identified) <- c("BarcodingR")

    return(output2.identified)


  #return(output.identified)

  }


  #bp(ref1,sampleSpeNames2,que1)


  #############################
  #Spp2<-as.factor(Spp)
  #rownames(ref1)<-Spp #sampleSpeNames
  #rownames(que1)<-queIDs



  fuzzyId<-function(ref1,Spp2,que1){
    #library(class)
    knn1<-knn(ref1, que1, cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    #knn1
    spe.Identified<-as.character(knn1)




    #rownames(ref1)<-sampleSpeNames

    # Ref<-ref1

    FMFtheta12<-function(Ref){

     # if (class(seqsRef)!="matrix")

      # stop("input should be an object of matrix!")


        ### 2.1 dealing with input DNA data!
        ### 2.2 calculate species center vectors
        morph.spe<-gsub(".+,","",rownames(Ref)) # remove sequence ID before ","
        #no.morph.spe<-length(unique(morph.spe))
        #species.centers<-aggregate(scale(Ref),by=list(morph.spe),FUN="mean")
        species.centers<-aggregate(Ref,by=list(morph.spe),FUN="mean")
        #class(species.centers)
        #head(species.centers)
        #dim(species.centers)


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

        #seqInOneSpe<-Ref[grep(levels(f)[1], Spp, value = FALSE,fixed = TRUE),]
        seqInOneSpe<-Ref[Spp %in% levels(f)[1],]
        #dim(seqInOneSpe)==NULL
        ifelse(popSize.PS[1]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean)) ###
        ifelse(popSize.PS[1]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean)) ###  centroid.spe0 -

        length.sites<-length(centroid.spe)
        ifelse(popSize.PS[1]==1,intra.dist<-0,intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe)) ### to the centroid
        ifelse(popSize.PS[1]==1,theta1.tmp[1]<-0,theta1.tmp[1]<-3*sd(intra.dist)) ### theta1 is sligthly different from
        #ifelse(popSize.PS[1]==1,theta1.tmp[1]<-0,theta1.tmp[1]<-max(intra.dist)) ### theta1 is sligthly different from

        ######



        for(i in 2:length(levels(f))){
          #i<-2
          #   i<-9

          #i=3
          #cat(paste("i=",i),"\n")
          # cat("\n")


          #seqInOneSpe<-sDNAbin[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),]
          #seqInOneSpe<-Ref[grep(levels(f)[i], Spp, value = FALSE,fixed = TRUE),]
          seqInOneSpe<-Ref[Spp %in% levels(f)[i],]
          ifelse(popSize.PS[i]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean))
          #ifelse(popSize.PS[i]==1,centroid.spe<-seqInOneSpe,centroid.spe<-apply(seqInOneSpe,2,mean))
          ifelse(popSize.PS[i]==1,centroid.spe<-c(centroid.spe,seqInOneSpe),centroid.spe<-c(centroid.spe,apply(seqInOneSpe,2,mean)))


          #length(seqInOneSpe)

          #ifelse(dim(seqInOneSpe)==NULL,theta1.tmp[i]<-0,theta1.tmp[i]<-max(dist(seqInOneSpe))) # error!
          #ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-max(dist(seqInOneSpe)))

          ifelse(popSize.PS[i]==1,intra.dist<-eucl.dist.two.vect(seqInOneSpe,centroid.spe0),intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0))


          #intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0) ### to the centroid
          #ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-max(intra.dist)) ### theta1 is sligthly different from
          ifelse(popSize.PS[i]==1,theta1.tmp[i]<-0,theta1.tmp[i]<-3*sd(intra.dist)) ### theta1 is sligthly different from


          #theta1.tmp[i]<-max(dist(seqInOneSpe))


        }




        centroid.spe.matrix<-t(array(centroid.spe,dim=c(length.sites,length(centroid.spe)%/%length.sites)))


        ### 1.2 calculate theta2 for each pairt of species


        ###################
        #codes<-out.somu$out.som.unique$codes  ### seqsRef$codes
        codes<-centroid.spe.matrix

        # dim(codes)


        theta12.all<-data.frame(list.spe=list.spe,PS=pairs[,1],NN=pairs[,2],stringsAsFactors=TRUE)

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
    FMF.theta12<-FMFtheta12(ref1)

    #mpattern<-".+,Noctuidae_"
    #mpattern<-"Noctuidae_"
    #mpattern<-"Noctuidae_"
    #FMF.theta12$list.spe<-gsub(mpattern,"",FMF.theta12$list.spe) # remove seqs names before "," (incl.",")
    #FMF.theta12$list.spe


    FMF.que<-numeric(dim(que1)[1])

    for(j in 1:dim(que1)[1]){
      #j<-81

    #seqInOneSpe<-ref1[grep(spe.Identified[j], FMF.theta12$list.spe, value = FALSE,fixed = TRUE),]
    seqInOneSpe<-ref1[grep(spe.Identified[j], rownames(ref1), value = FALSE,fixed = TRUE),]
    #class(FMF.theta12$list.spe)
    k<-match(x=spe.Identified[j], table=FMF.theta12$list.spe)

    #ifelse(FMF.theta12$popSize.PS[k]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean))
    ifelse(FMF.theta12$popSize.PS[k]==1,centroid.spe0<-seqInOneSpe,centroid.spe0<-apply(seqInOneSpe,2,mean))


    que2PS.dist<-eucl.dist.two.vect(que1[j,],centroid.spe0)
    #que2PS.dist<-0
    #intra.dist<-apply(seqInOneSpe,1,eucl.dist.two.vect,v2=centroid.spe0)

    xtheta12<-c(que2PS.dist,FMF.theta12$theta1[k],FMF.theta12$theta2[k])
    FMF.que[j]<-FMF(xtheta12)
    }

    #FMF.theta12<-FMFtheta12(ref1)

    #sub.FMF.theta12<-subset(FMF.theta12,popSize.PS>1)
    #average.theta1<-mean(sub.FMF.theta12$theta1)

    #### calculate success.rates.ref

    knn1<-knn(ref1, ref1, cl=Spp2, k = 1, l = 0, prob = FALSE, use.all = TRUE)
    #knn1
    spe.Identified.ref<-as.character(knn1)



    spe.morph.Identified<-data.frame(Spp,spe.Identified.ref,stringsAsFactors=TRUE)

    matches<-apply(spe.morph.Identified,2,strings.equal,str2=spe.morph.Identified[,2])


    matches<-colSums(matches,dims = 1)

    success.rates.ref<-matches[1]/matches[2]
    names(success.rates.ref)<-NULL



    output.identified<-data.frame(queIDs=queIDs,
                                  spe.Identified=spe.Identified,
                                  FMF=FMF.que,stringsAsFactors=TRUE)

    out<-list(success.rates.ref=success.rates.ref,
              #output_identified.ref=output.identified.ref,
              output_identified=output.identified)


    #Bayesian.prob=Bayesian.prob)

    # return(output.identified)

    class(out) <- c("BarcodingR")

    return(out)

    ##################
   # output.identified<-data.frame(queIDs=queIDs,
     #                             spe.Identified=spe.Identified,
     #                             FMF=FMF.que)
    #Bayesian.prob=Bayesian.prob)

   # return(output.identified)




  }

  #fuzzy.Id<-fuzzyId(ref1,Spp2,que1)


  Bayesian<-function(ref,Spp,que){

    rq<-rbind(ref,que)
    rq<-rq[,seg.sites(rq)] # remove constant sites!

    queIDs<-attr(que,"dimnames")[[1]]


    ref.constant.sites.removed<-rq[1:dim(ref)[1],]
    que.constant.sites.removed<-rq[-(1:dim(ref)[1]),]

    #rownames(ref.constant.sites.removed)
    #rownames(que.constant.sites.removed)

    ref2<-as.character(ref.constant.sites.removed)
    ref2<-as.data.frame(ref2,stringsAsFactors=TRUE)
    que2<-as.character(que.constant.sites.removed)
    que2<-as.data.frame(que2,stringsAsFactors=TRUE)
    rq2<-rbind(ref2,que2)
    Spp<-c(Spp,queIDs)
    rq2$species<-as.factor(Spp)
    ref3<-rq2[1:dim(ref2)[1],]

    #rq2$species



    ref3$species<-with(ref3,as.factor(gsub(".+,","",sampleSpeNames)))

    #ref3$species

    que3<-rq2[-(1:dim(ref2)[1]),]
    del<-dim(que3)[2]
    que3<-que3[,-del]


    head(ref3);class(ref3);dim(ref3)

    rownames(ref3)<-1:dim(ref3)[1]

    Bayesian.trained <- naiveBayes(species ~ ., data = ref3)
    #Bayesian.trained <- naiveBayes(species ~ ., data = ref3,type = "raw")
    #model <- naiveBayes(species ~ ., data = training.set)
    #attributes(Bayesian.trained)
    #attr(Bayesian.trained,"apriori")[[1]]

    spe.inferred<-predict(Bayesian.trained, que3)
    spe.inferred.prob<-predict(Bayesian.trained, que3, type = "raw")
    #dim(que3)
    #attributes(spe.inferred)
    Bayesian.prob<-apply(spe.inferred.prob,1,max)

    #spe.inferred<-predict(Bayesian.trained, que3,type = "raw")
    #spe.inferred
    spe.inferred<-as.character(spe.inferred)



    output.identified<-data.frame(queIDs=queIDs,
                                  spe.Identified=spe.inferred,
                                  Bayesian.prob=Bayesian.prob,stringsAsFactors=TRUE)

    #class(output.identified) <- c("BarcodingR")


    out<-list(output_identified=output.identified)


    #Bayesian.prob=Bayesian.prob)

    # return(output.identified)

    class(out) <- c("BarcodingR")



    return(out)
  }

  #Bayesian(ref,Spp,que)

 # bp(ref1,sampleSpeNames2,ref1)

  #bpNewTraining, bpUseTrained
  method <- pmatch(method, c("bpNewTraining", ### train and identify at the same time
                             "bpNewTrainingOnly", ### just train the model for later use
                             "bpUseTrained",
                             "fuzzyId",
                             "Bayesian",
                             "all"))
  if (is.na(method))
    stop("invalid method")
  if (method == -1)
    stop("ambiguous method")

  if (method == 1){  ### "bpNewTraining", ### train and identify at the same time
    sampleSpeNames<-attr(ref,"dimnames")[[1]]
    queIDs<-attr(que,"dimnames")[[1]]

    mpattern<-".+,"
    #mpattern<-".+,Noctuidae_"
    #mpattern<-"Noctuidae_"
    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    Spp
    #length(unique(Spp))
    sampleSpeNames2<-class.ind(Spp)

    ### rq
    rq<-rbind(ref,que)
    rq<-rq[,seg.sites(rq)] # remove constant sites!

    ref1<-rq[1:dim(ref)[1],]
    que1<-rq[-(1:dim(ref)[1]),]

    ref1<-digitize.DNA(ref1)
    #names(ref)<-sampleSpeNames
    que1<-digitize.DNA(que1)

    out<-bp(ref1,sampleSpeNames2,method = "bpNewTraining",que1)

  }


  if (method == 2){  ###"bpNewTrainingOnly" ### just train the model for later use
    sampleSpeNames<-attr(ref,"dimnames")[[1]]
    #queIDs<-attr(que,"dimnames")[[1]]

    mpattern<-".+,"

    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    Spp

    sampleSpeNames2<-class.ind(Spp)

    ref1<-digitize.DNA(ref)
    #que1<-digitize.DNA(que1)
    out<-bp(ref1,sampleSpeNames2,method = "bpNewTrainingOnly",que1)



  }



  if (method == 3){ ###"bpUseTrained"
    #que <-fasta2DNAbin("que.fas")
    que1<-digitize.DNA(que)

    out<-bp(ref1,sampleSpeNames2,method = "bpUseTrained",que1)
  }
  if (method == 4){###  "fuzzyId"


    rq<-rbind(ref,que)
    rq<-rq[,seg.sites(rq)] # remove constant sites!

    ref1<-rq[1:dim(ref)[1],]
    que1<-rq[-(1:dim(ref)[1]),]

    ref1<-digitize.DNA(ref1)
    #names(ref)<-sampleSpeNames
    que1<-digitize.DNA(que1)


    Spp2<-as.factor(Spp)
    rownames(ref1)<-Spp #sampleSpeNames
    queIDs<-attr(que,"dimnames")[[1]]
    rownames(que1)<-queIDs



    out<-fuzzyId(ref1,Spp2,que1)
  }



  if (method == 5)
    out<-Bayesian(ref,Spp,que)

  return(out)

}




#setwd("C:/R/myRprojects/SpeDelimitation/SOFM")
#source("run.time.R")
#start.time<-Sys.time()

#####

#bsi<-barcoding.spe.identify(ref, que, method = "bpNewTraining")

#####




#setwd("C:/R/myRprojects/SpeDelimitation/SOFM")
#time.elapsed<-run.time(start.time)
#cat("Time used : H:M:S")
#time.elapsed

























