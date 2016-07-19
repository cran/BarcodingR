#' Sample Random Datasets from References (DNAbin)
#' 
#' @description Randomly sample reference data at different levels of taxon.
#' @param  ref Object of class "DNAbin" used as a reference dataset, which contains taxon information.
#' @param  sample.porp a numeric value between 0 and 1, indicating proportion of samples to draw 
#' at each level of taxon.
#' @param  sample.level a character string choosing from c("full","family","genus","species"). 
#' @return a list containing the selected samples and the samples left, in DNAbin format stored in a matrix or a list.
#' @keywords sample.ref
#' @export 
#' @note the ref must contain information on taxonomy, in format like, ">LS0909030M,Noctuidae_Himalaea_unica",
#' i.e., "seqID,family_genus_species", or ">LS0909030M,Himalaea_unica"; in case there is only one sample/individual
#' for a taxon level, this sample will be retained in ref.selected.
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA.
#' @references zhangab2008(at)mail.cnu.edu.cn;

#' @examples
#' 
#' data(TibetanMoth) 
#' data(pineMothITS2)
#' ref<-TibetanMoth
#' ref2<-pineMothITS2
#' out<-sample.ref(ref,sample.porp=0.5,sample.level="full")
#' out
#' out2<-sample.ref(ref2,sample.porp=0.5,sample.level="full")
#' out2
#' 
#' 

sample.ref<-function(ref,sample.porp=0.5,sample.level="full"){
  ###sample.level<-c("full","family","genus","species")
  ### functions used:
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
  if (class(ref)!="DNAbin")
    stop("seqs should be in DNAbin format!")
  
  ### 1. basic statistics for taxon information
  sampleSpeNames<-NAMES(ref)
  id<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="id")
  family<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="family") ### just return a list!
  family<-as.character(family)
  
  genus<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="genus")
  genus<-as.character(genus)
  
  species<-sapply(sampleSpeNames,taxonInfoExtraction,returnValue="species")
  species<-as.character(species)
  
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
                  species.list=species.list)
  
  ### 2. make random samples specified by the parameter "sample.level"
  ### 2.1 case "full"
  ref.selected<-NULL
  ref.left<-NULL

  if(sample.level=="full"){
    if(mode(ref)=="raw"){### case COI 
      sampled<-sample(dim(ref)[1], size=dim(ref)[1]*sample.porp)
      ref.selected<-ref[sampled,]
      ref.left<-ref[-sampled,]
      
      #ref.selected<-ref[sample(dim(ref)[1], size=dim(ref)[1]*sample.porp),]
      #ref.left<-ref[-sample(dim(ref)[1], size=dim(ref)[1]*sample.porp),]
      
    }else{ ### mode(ref)=="list" case ITS
      #ref<-ref2
      sampled<-sample(length(ref), size=length(ref)*sample.porp)
      ref.selected<-ref[sampled]
      ref.left<-ref[-sampled]
      
      #ref.selected<-ref[sample(length(ref), size=length(ref)*sample.porp)]
      #ref.left<-ref[-sample(length(ref), size=length(ref)*sample.porp)]
      #names(ref)
    }
    
  }
  
  ### 2.2 case "family"
   if(sample.level=="family"){
     if(family[1]=="NA"){ ### LS0909030M,Himalaea_unica ## no family information!
       stop("Your data have no family information!!!")
     }else{
    sampleSpeNames<-NAMES(ref)
    mpattern<-".+,"
    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
     
    ### loop to deal with all families:
    ### initinize ref.selected and ref.left using the first loop
    ref.selected<-NULL
    ref.left<-NULL
    
    if(mode(ref)=="raw"){### case COI 
      taxon.name<-paste(taxonStat$family.list[1], "_", sep = "")
      each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),] 
      
      #each.taxon<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
    if(dim(each.taxon)[1]>1){
    sampled<-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp)
    ref.selected<-each.taxon[sampled,]
    ref.left<-each.taxon[-sampled,]
    }else{
      ref.selected<-each.taxon
    }
    #ref.selected<-each.taxon[sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
    #ref.left<-each.taxon[-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
    
    
    if(length(taxonStat$family.list)>1){
    for(i in 2:length(taxonStat$family.list)){
    ### i<-2
      cat("i:",i,"\n")
    ### 2.2.1 extract all sequences for each family, and choose specified proportion of that
    ### each family into ref.selected, the left into ref.left
    #seqInOneSpe<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
      
      taxon.name<-paste(taxonStat$family.list[i], "_", sep = "")
      each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),]
      
      #each.taxon<-ref[grep(taxonStat$family.list[i], Spp, value = FALSE,fixed = TRUE),] 
      
      #cat("no.sample.each taxon:",dim(each.taxon)[1],"\n")
      
      
      
      if(dim(each.taxon)[1]>1){
      sampled<-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp)
      ref.selected.tmp<-each.taxon[sampled,]
      ref.left.tmp<-each.taxon[-sampled,]
      
      #cat("no.sample.ref.selected.tmp:",dim(ref.selected.tmp)[1],"\n")
      
      #cat("no.sample.ref.left.tmp:",dim(ref.left.tmp)[1],"\n")
      
      #ref.selected.tmp<-each.taxon[sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
      #ref.left.tmp<-each.taxon[-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
    
      ref.selected<-rbind(ref.selected,ref.selected.tmp)
      #ref.left<-rbind(ref.left,ref.left.tmp)
      if(length(ref.left)>0&length(ref.left.tmp)>0){
        ref.left<-rbind(ref.left,ref.left.tmp)}else{ref.left<-ref.left.tmp }
      }else{
        ref.selected.tmp<-each.taxon
        ref.selected<-rbind(ref.selected,ref.selected.tmp)
        
       # cat("just one sample!\n")
        
      }
      
    }### the end of i for-loop
    }### the end of if-length(taxonStat$family.list) loop
    }else{ ### mode(ref)=="list" case ITS
      taxon.name<-paste(taxonStat$family.list[1], "_", sep = "")
      each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),]
      #each.taxon<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE)] 
      if(length(each.taxon)>1){
      sampled<-sample(length(each.taxon), size=length(each.taxon)*sample.porp)
      ref.selected<-each.taxon[sampled]
      ref.left<-each.taxon[-sampled]
      }else{
        ref.selected<-each.taxon 
      }
      #ref.selected<-each.taxon[sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
      #ref.left<-each.taxon[-sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
      
      if(length(taxonStat$family.list)>1){
        for(i in 2:length(taxonStat$family.list)){
          ### i<-2
          
          ### 2.2.1 extract all sequences for each family, and choose specified proportion of that
          ### each family into ref.selected, the left into ref.left
          #seqInOneSpe<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
          taxon.name<-paste(taxonStat$family.list[i], "_", sep = "")
          each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),]
          #each.taxon<-ref[grep(taxonStat$family.list[i], Spp, value = FALSE,fixed = TRUE)]
          if(length(each.taxon)>1){
          sampled<-sample(length(each.taxon), size=length(each.taxon)*sample.porp)
          ref.selected.tmp<-each.taxon[sampled]
          ref.left.tmp<-each.taxon[-sampled]
          
          
          #ref.selected.tmp<-each.taxon[sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
          #ref.left.tmp<-each.taxon[-sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
          
          
          ref.selected<-c(ref.selected,ref.selected.tmp)
          #ref.left<-c(ref.left,ref.left.tmp)
          
          
          if(length(ref.left)>0&length(ref.left.tmp)>0){
            ref.left<-c(ref.left,ref.left.tmp)}else{ref.left<-ref.left.tmp}
          
          
          
          }else{
            ref.selected.tmp<-each.taxon
            ref.selected<-c(ref.selected,ref.selected.tmp)
            
          }
          #ref.selected<-rbind(ref.selected,ref.selected.tmp)
          #ref.left<-rbind(ref.left,ref.left.tmp)
          
        }### the end of i for-loop
      }### the end of if-length(taxonStat$family.list) loop
      
      
    }
    
    
  }
   }
  

### 2.3 case "genus"
  if(sample.level=="genus"){
    #taxonStat$family.list[1]
    #class(taxonStat$family.list)
    #attributes(taxonStat)
    #taxonStat$genera.list[1]
    sampleSpeNames<-NAMES(ref)
    mpattern<-".+,"
    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    
    ### loop to deal with all genera:
    ### initinize ref.selected and ref.left using the first loop
    ref.selected<-NULL
    ref.left<-NULL
    
    if(mode(ref)=="raw"){### case COI 
      #each.taxon<-ref[grep(taxonStat$genera.list[1], Spp, value = FALSE,fixed = TRUE),] 
      taxon.name<-paste(taxonStat$genera.list[1], "_", sep = "")
      each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),] 
      
      if(dim(each.taxon)[1]>1){
      
      sampled<-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp)
      ref.selected<-each.taxon[sampled,]
      ref.left<-each.taxon[-sampled,]
      }else{
        ref.selected<-each.taxon
      }
      
      
      #ref.selected<-each.taxon[sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
      #ref.left<-each.taxon[-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
      
      
      if(length(taxonStat$genera.list)>1){
        for(i in 2:length(taxonStat$genera.list)){
          ### i<-44
          cat("i:",i,"\n")
          ### 2.2.1 extract all sequences for each family, and choose specified proportion of that
          ### each family into ref.selected, the left into ref.left
          #seqInOneSpe<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
          
          taxon.name<-paste(taxonStat$genera.list[i], "_", sep = "")
          each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),] 
          #each.taxon<-ref[grep(taxonStat$genera.list[i], Spp, value = FALSE,fixed = TRUE),] 
          
          if(dim(each.taxon)[1]>1){
          sampled<-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp)
          
          #ref.selected.tmp<-each.taxon[sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
          #ref.left.tmp<-each.taxon[-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
          
          ref.selected.tmp<-each.taxon[sampled,]
          ref.left.tmp<-each.taxon[-sampled,]
          
          ref.selected<-rbind(ref.selected,ref.selected.tmp)
          #ref.left<-rbind(ref.left,ref.left.tmp)
          if(length(ref.left)>0&length(ref.left.tmp)>0){
            ref.left<-rbind(ref.left,ref.left.tmp)}else{ref.left<-ref.left.tmp }
          
          }else{
            ref.selected.tmp<-each.taxon
            ref.selected<-rbind(ref.selected,ref.selected.tmp)
          }
          #ref.selected<-rbind(ref.selected,ref.selected.tmp)
          #ref.left<-rbind(ref.left,ref.left.tmp)
          
          
        }### the end of i for-loop
      }### the end of if-length(taxonStat$family.list) loop
    }else{ ### mode(ref)=="list" case ITS
      #taxon.name<-paste(taxonStat$genera.list[i], "_", sep = "")
      taxon.name<-paste(taxonStat$genera.list[1], "_", sep = "")
      #each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),] 
      each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE)] 
      
      #each.taxon<-ref[grep(taxonStat$genera.list[1], Spp, value = FALSE,fixed = TRUE)] 
      if(length(each.taxon)>1){
      sampled<-sample(length(each.taxon), size=length(each.taxon)*sample.porp)
      ref.selected<-each.taxon[sampled]
      ref.left<-each.taxon[-sampled]
      }else{
        ref.selected<-each.taxon
      }
      
      #ref.selected<-each.taxon[sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
      #ref.left<-each.taxon[-sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
      
      
      if(length(taxonStat$genera.list)>1){
        for(i in 2:length(taxonStat$genera.list)){
          ### i<-2
          
          ### 2.2.1 extract all sequences for each family, and choose specified proportion of that
          ### each family into ref.selected, the left into ref.left
          #seqInOneSpe<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
          taxon.name<-paste(taxonStat$genera.list[i], "_", sep = "")
          #each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE),] 
          each.taxon<-ref[grep(taxon.name, Spp, value = FALSE,fixed = TRUE)]
          #each.taxon<-ref[grep(taxonStat$genera.list[i], Spp, value = FALSE,fixed = TRUE)] 
          if(length(each.taxon)>1){
          sampled<-sample(length(each.taxon), size=length(each.taxon)*sample.porp)
          ref.selected.tmp<-each.taxon[sampled]
          ref.left.tmp<-each.taxon[-sampled]
          
          ref.selected<-c(ref.selected,ref.selected.tmp)
          #ref.left<-c(ref.left,ref.left.tmp)
          
          
          if(length(ref.left)>0&length(ref.left.tmp)>0){
            ref.left<-c(ref.left,ref.left.tmp)}else{ref.left<-ref.left.tmp}
          
          
          
          }else{
            ref.selected.tmp<-each.taxon
            ref.selected<-c(ref.selected,ref.selected.tmp)
            
          }
          #ref.selected.tmp<-each.taxon[sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
          #ref.left.tmp<-each.taxon[-sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
          
          
          #ref.selected<-rbind(ref.selected,ref.selected.tmp)
          #ref.left<-rbind(ref.left,ref.left.tmp)
          
          
        }### the end of i for-loop
      }### the end of if-length(taxonStat$family.list) loop
      
      
    }
    
    
  }

  #}
  ### 2.4 case "species"
  if(sample.level=="species"){
    sampleSpeNames<-NAMES(ref)
    mpattern<-".+,"
    Spp<-gsub(mpattern,"",sampleSpeNames) # remove seqs names before "," (incl.",")
    
    Spp_tmp<-paste(Spp,"_",sep = "")
    
    ### loop to deal with all species:
    ### initinize ref.selected and ref.left using the first loop
    ref.selected<-NULL
    ref.left<-NULL
    
    if(mode(ref)=="raw"){### case COI 
      taxon.name<-paste(taxonStat$species.list[1], "_", sep = "")
      each.taxon<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE),]
      
      #each.taxon<-ref[grep(taxonStat$species.list[1], Spp, value = FALSE,fixed = TRUE),] 
      #cat("no.sample.each taxon:",dim(each.taxon)[1],"\n")
      
      if(dim(each.taxon)[1]>1){
      
      sampled<-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp)
      ref.selected<-each.taxon[sampled,]
      ref.left<-each.taxon[-sampled,]
      
      #cat("no.sample.ref.selected.tmp:",dim(ref.selected.tmp)[1],"\n")
      
      #cat("no.sample.ref.left.tmp:",dim(ref.left.tmp)[1],"\n")
      
      }else{
        ref.selected<-each.taxon
        
        cat("just one sample!\n")
      }
      #ref.selected<-each.taxon[sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
      #ref.left<-each.taxon[-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
      
      
      if(length(taxonStat$species.list)>1){
        for(i in 2:length(taxonStat$species.list)){
          ### i<-3
          cat("i:",i,"\n")
          ### 2.2.1 extract all sequences for each family, and choose specified proportion of that
          ### each family into ref.selected, the left into ref.left
          #seqInOneSpe<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
          
          taxon.name<-paste(taxonStat$species.list[i], "_", sep = "")
          each.taxon<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE),]
          
          #each.taxon<-ref[grep(taxonStat$species.list[i], Spp, value = FALSE,fixed = TRUE),] 
          #cat("no.sample.each taxon:",dim(each.taxon)[1],"\n")
          
          if(dim(each.taxon)[1]>1){
          #dim(each.taxon)
          sampled<-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp)
          ref.selected.tmp<-each.taxon[sampled,]
          ref.left.tmp<-each.taxon[-sampled,]
          
          #cat("no.sample.ref.selected.tmp:",dim(ref.selected.tmp)[1],"\n")
          #cat("no.sample.ref.left.tmp:",dim(ref.left.tmp)[1],"\n")
          
          #ref.selected.tmp<-each.taxon[sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
          #ref.left.tmp<-each.taxon[-sample(dim(each.taxon)[1], size=dim(each.taxon)[1]*sample.porp),]
          #class()
          ref.selected<-rbind(ref.selected,ref.selected.tmp)
          if(length(ref.left)>0&length(ref.left.tmp)>0){
          ref.left<-rbind(ref.left,ref.left.tmp)}else{ref.left<-ref.left.tmp }
          }else{
            ref.selected.tmp<-each.taxon
            ref.selected<-rbind(ref.selected,ref.selected.tmp)
            
            #cat("just one sample!\n")
          }
          
          
          #cat("ref.selected:",dim(ref.selected)[1],"\n")
          #cat("ref.left:",dim(ref.left)[1],"\n")
          
        }### the end of i for-loop
        
        
      }### the end of if-length(taxonStat$family.list) loop
    }else{ ### mode(ref)=="list" case ITS
      taxon.name<-paste(taxonStat$species.list[1], "_", sep = "")
      #each.taxon<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE),]
      each.taxon<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE)]
      
      #each.taxon<-ref[grep(taxonStat$species.list[1], Spp, value = FALSE,fixed = TRUE)] 
      if(length(each.taxon)>1){
      sampled<-sample(length(each.taxon), size=length(each.taxon)*sample.porp)
      ref.selected<-each.taxon[sampled]
      ref.left<-each.taxon[-sampled]
      }else{
        ref.selected<-each.taxon
      }
      #ref.selected<-each.taxon[sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
      #ref.left<-each.taxon[-sample(length(each.taxon), size=length(each.taxon)*sample.porp)]
      
      
      
      if(length(taxonStat$species.list)>1){
        for(i in 2:length(taxonStat$species.list)){
          ### i<-2
          
          ### 2.2.1 extract all sequences for each family, and choose specified proportion of that
          ### each species into ref.selected, the left into ref.left
          #seqInOneSpe<-ref[grep(taxonStat$family.list[1], Spp, value = FALSE,fixed = TRUE),] 
          taxon.name<-paste(taxonStat$species.list[i], "_", sep = "")
          #each.taxon<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE),]
          each.taxon<-ref[grep(taxon.name, Spp_tmp, value = FALSE,fixed = TRUE)]
          
          #each.taxon<-ref[grep(taxonStat$species.list[i], Spp, value = FALSE,fixed = TRUE)] 
          if(length(each.taxon)>1){
          sampled<-sample(length(each.taxon), size=length(each.taxon)*sample.porp)
          
          ref.selected.tmp<-each.taxon[sampled]
          ref.left.tmp<-each.taxon[-sampled]
          
          #ref.selected<-list(ref.selected,ref.selected.tmp)
          #ref.left<-list(ref.left,ref.left.tmp)
          
          #ref.selected<-rbind(ref.selected,ref.selected.tmp)
          #ref.left<-rbind(ref.left,ref.left.tmp)
          ref.selected<-c(ref.selected,ref.selected.tmp)
          #ref.left<-c(ref.left,ref.left.tmp)
          
          if(length(ref.left)>0&length(ref.left.tmp)>0){
            ref.left<-c(ref.left,ref.left.tmp)}else{ref.left<-ref.left.tmp}
          
          
          }else{
            ref.selected.tmp<-each.taxon
            ref.selected<-c(ref.selected,ref.selected.tmp)
          }
        }### the end of i for-loop
      }### the end of if-length(taxonStat$family.list) loop
      
      
    }
    
    
  }

  #}
  
  out<-list(ref.selected=ref.selected,ref.left=ref.left)
  return(out)
  
}


