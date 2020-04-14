#' Save Identifications
#' 
#' @description Output identified results to an outfile in temporty directory (found by tempdir() function).
#' 
#' @param  outfile character string to indicate outfile name.
#' @param  ids object of class "BarcodingR", which contains identified taxon information. 
#' 
#' @return no value returned,but an output file.
#' @keywords save.ids
#' @export 
#' @author Ai-bing ZHANG, PhD. CNU, Beijing, CHINA.
#' @references zhangab2008(at)mail.cnu.edu.cn
#' @seealso barcoding.spe.identify()
#' @examples
#' 
#' 
#' data(TibetanMoth) 
#' ref<-as.DNAbin(as.character(TibetanMoth[1:50,]))
#' que<-as.DNAbin(as.character(TibetanMoth[50:60,]))
#' bsi<-barcoding.spe.identify(ref, que, method = "fuzzyId")
#' bsi
#' save.ids(outfile="identified.txt",bsi)

save.ids<-function(outfile="identified.txt",ids){
  if(class(ids)!="BarcodingR")
    stop("A BarcodingR object is required for ids!!!")
  
  #attributes(bsi)
  
  Rhome<-tempdir() ### 2020/4/13 21:06:46
  fileName<-outfile
  fileName<-paste(Rhome,fileName,sep = "/")### 2020/4/13 21:06:46
  #fileName<-paste(Rhome,fileName,sep = "\\")### 2020/4/13 21:06:46
  fileName
  #fileName<-paste("simulation",i,sep = "")
  #fileName<-paste(fileName,".RData",sep = "")
  #fileName
  outfile<-fileName
  x<-ids$output_identified
  #x<-bsi$output_identified
  write.table(x, file = outfile, append = FALSE, quote = FALSE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
}

#save.ids(outfile="out.txt",bsi)
