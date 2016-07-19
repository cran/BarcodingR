#' Summarize Reference Data
#' 
#' @description Summarize taxon information, sequence statistics,barcodes numbers per species for reference dataset. 
#' 
#' @param  ref object of class "DNAbin" used as a reference dataset, which contains taxon information
#' @param  taxonStat logic value to indicate whether the item is calculated
#' @param  seqStat logic value to indicate whether the item is calculated
#' @param  barcodeStat logic value to indicate whether the item is calculated
#' @return a list containing taxon statistics, sequence statistics, population parameters,barcoding statistics ()
#' @keywords summarize ref
#' @export 
#' @author Ai-bing ZHANG, Meng-di HAO, CNU, Beijing, CHINA.
#' @references zhangab2008(at)mail.cnu.edu.cn./zhangab2008(at)gmail.com.
#' @examples
#' 
#' data(TibetanMoth) 
#' s.r<-summarize.ref(TibetanMoth,taxonStat=TRUE,seqStat=TRUE,barcodeStat=TRUE)
#' s.r

summarize.ref<-function(ref,taxonStat=TRUE,seqStat=TRUE,barcodeStat=TRUE){

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
  
  if (class(ref)!="DNAbin")
    stop("seqs should be in DNAbin format!")
  
  ### 1. taxon statistics:
  if (taxonStat==TRUE){
    
    
    #sampleSpeNames<-attr(ref,"dimnames")[[1]]
    sampleSpeNames<-NAMES(ref)
    #ifelse(length(sampleSpeNames)!=0,sampleSpeNames<-sampleSpeNames,sampleSpeNames<-names(ref))
    

    taxonInfoExtraction<-function(seqLables,returnValue="id"){
      
      if(class(seqLables)!="character")
        stop("seqLables is not character!")
      
      id<-strsplit(seqLables, ",")[[1]][1]
      taxon<-strsplit(seqLables, ",")[[1]][2]
      taxa<-strsplit(taxon, "_")
      taxa<-taxa[[1]]
      
      
      if(length(taxa)==3){ ### with family information, "LS0909030M,Noctuidae_Himalaea_unica"
        
        family<-taxa[1]
        genus<-paste(taxa[1],taxa[2],sep="_")
        species<-paste(taxa[2],taxa[3],sep="_")
        
        #ti<-c(id,family,genus,species)
        
      }else{ ### no family information, "LS0909030M,Himalaea_unica"
        family<-"NA"
        genus<-taxa[1]
        #genus<-paste(taxa[1],taxa[2],sep="_")
        species<-paste(taxa[1],taxa[2],sep="_")
        #ti<-c(id,genus,species)
        
      }
      #ti<-data.frame(IDs=id,families=family,genera=genus,species=species)
      #ti<-c(id,family,genus,species)
      #return(ti)
      #family<-as.character(family)
      if(returnValue=="id") return(id)
      if(returnValue=="family") return(family)
      if(returnValue=="genus") return(genus)
      if(returnValue=="species") return(species)
      
    }
    
    #id<-apply(sampleSpeNames,1,FUN=taxonInfoExtraction,returnValue="id")
    #family<-apply(sampleSpeNames,1,FUN=taxonInfoExtraction,returnValue="family")
    id<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="id")
    family<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="family") ### just return a list!
    family<-as.character(family)
    genus<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="genus")
    genus<-as.character(genus)
    
    species<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="species")
    species<-as.character(species)
    
    sample.sizes<-list(family=table(family),genus=table(genus),species=table(species))
    
    
    if(family[1]!="NA"){ ### LS0909030M,Noctuidae_Himalaea_unica
      IDs<-id
      no.samples<-length(IDs)
      
      family.list<-unique(family)
      no.family<-length(unique(family))
      
      genera.list<-unique(genus)
      no.genera<-length(unique(genus))
      
      species.list<-unique(species)
      no.species<-length(unique(species))
      
      basic.stat<-c(no.samples,no.family,no.genera,no.species)
      names(basic.stat)<-c("no.samples","no.family","no.genera","no.species")
      
      taxonStat<-list(basic.stat=basic.stat,
                      family.list=family.list,
                      genera.list=genera.list,
                      species.list=species.list,
                      sample.sizes=sample.sizes
                      )

    }else{### LS0909030M,Himalaea_unica
      IDs<-id
      no.samples<-length(IDs)
      
      #family.list<-unique(taxa[2,])
      #no.family<-length(unique(taxa[2,]))
      
      genera.list<-unique(genus)
      no.genera<-length(unique(genus))
      
      species.list<-unique(species)
      no.species<-length(unique(species))
      
      basic.stat<-c(no.samples,no.genera,no.species)
      names(basic.stat)<-c("no.samples","no.genera","no.species")
      
      taxonStat<-list(basic.stat=basic.stat,
                      genera.list=genera.list,
                      species.list=species.list,
                      sample.sizes=sample.sizes
                )
      
    } ### the end of else
  }else{taxonStat<-NULL}
  
  if (seqStat==TRUE){
    no.seqs<-dim(ref)[1]
    length.seqs<-dim(ref)[2]
    ifelse(length(no.seqs)!=0,no.seqs<-no.seqs,no.seqs<-length(ref))
    
    ifelse(length(length.seqs)!=0,length.seqs<-length.seqs,length.seqs<-mean(as.numeric(summary(ref)[,1])))
      seqStat<-c(no.seqs,length.seqs)
    names(seqStat)<-c("no.seqs","length.seqs")
  }else{seqStat<-NULL}
  
  if (barcodeStat==TRUE){
    #sampleSpeNames<-attr(ref,"dimnames")[[1]]
    sampleSpeNames<-NAMES(ref)
       mpattern<-".+,"
      Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    #Spp
      Spp_tmp<-paste(Spp,"_",sep = "")
      #f<-factor(Spp)
      f<-factor(Spp_tmp) 
    pop.size<-numeric(length(levels(f)))
    
    for(i in 1:length(levels(f))){
      #i<-1
      taxon.name<-paste(levels(f)[i], "_", sep = "")
      #samplesInOneSpe<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE),]
      #samplesInOneSpe<-list(grep(levels(f)[i], Spp, value = TRUE))
      samplesInOneSpe<-list(grep(levels(f)[i], Spp_tmp, value = TRUE))
      
      
      #DNAsamplesInOneSpe<-minput.fas[grep(levels(f)[i], Spp, value = F),]#updata2014-4-24
        pop.size[i]<-length(samplesInOneSpe[[1]])
      
    }
    
    min.popSize<-min(pop.size)
    max.popSize<-max(pop.size)
    mean.popSize<-mean(pop.size)
     basic.stat<-c(min.popSize,max.popSize,mean.popSize)
    names(basic.stat)<-c("min.popSize","max.popSize","mean.popSize")
    
    
    barcodeStat<-list(basic.stat=basic.stat,
                      pop.size=data.frame(species=levels(f),pop.size=pop.size)
                )
    }else{barcodeStat<-NULL}
  
  
  ### 2. sequence statistics:
  ### 3. population parameters estimated:
  ### 4. barcode statistics: (interspecific, intraspecific distance, draw barcoding gaps)
  ###    barcodes/species (min,max,mean/spe)
  
  out<-list(taxonStat=taxonStat,
            seqStat=seqStat,
            barcodeStat=barcodeStat
    )
  return(out)
} 




















