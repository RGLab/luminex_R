# show method:
#   Shows class, slots, number of analytes, total number of measures
#
setMethod("show", "BAMAObject", function(object){
  cat("An object of class BAMAObject with",nrow(assayData(object)),"analytes:","\n")
  cat("\t", as.character(head(assayData(object)$name, 3)),"...", as.character(tail(assayData(object)$name, 3)),"\n")
  cat(length(unlist(exprs(object), use.names=FALSE)), "measures of expression in", nrow(pData(object)),"wells.","\n")
  cat("And slots:", names(getSlots("BAMAObject")),"\n")
})

#TODO: Add setters

# pData method:
#   phenoData accessor
setMethod("pData", "BAMAObject", function(object){
  return(object@phenoData)
})

# assayData
#   assayData accessor
setGeneric("assayData", function(object) standardGeneric("assayData"))
setMethod("assayData", "BAMAObject", function(object){
  return(object@assayData)
})

# summary
#   assayData accessor
setMethod("summary", "BAMAObject", function(object){
  return(object@summary)
})

# exprs method:
#   exprs accessor
setMethod("exprs", "BAMAObject", function(object){
  return(object@exprs)
})

# getMFI method
#   Calculates the MFI for a BAMAObject
#   OUTPUT: A matrix with analytes as rows and wells as columns
setGeneric("getMFI", function(object) standardGeneric("getMFI"))
setMethod("getMFI", "BAMAObject", function(object){
  matMFI<-sapply(exprs(object), function(x){sapply(x, median)})
  return(matMFI)
})

# subset method
#   subsets the object for one sample
subset.BAMAObject<-function(object, sample){
  pD<-pData(object)[which(pData(object)$ID==sample),]
  aD<-assayData(object)
  exprs<-list(exprs(object)[[sample]])
  names(exprs)<-sample #so the exprs retains it's named list of lists structure
  nBama<-new("BAMAObject", phenoData=pD, assayData=aD, exprs=exprs)
}

# as.data.frame method
#   Gather all available information into a big data.frame
setAs("BAMAObject", "data.frame", function(from){
  #all info is actually available from @exprs
  nSample<-nrow(pData(from))
  nAnalytes<-nrow(assayData(from))
  sampleVec<-c()
  BIDVec<-c()
  RP1Vec<-c()
  for(sample in 1:nSample)
  {
    for(analyte in 1:nAnalytes)
    {
      fluo<-exprs(from)[[sample]][[analyte]]
      RP1Vec<-c(RP1Vec,fluo)
      BIDVec<-c(BIDVec, rep(assayData(from)[analyte,"ID"], length(fluo)))
      sampleVec<-c(sampleVec, rep(as.character(pData(from)[sample,"ID"]), length(fluo)))
    }
  #TODO: lapply (ies?)
  }
  to<-data.frame(BID=BIDVec, RP1=RP1Vec, Well=sampleVec)
  return(to)
})

as.data.frame.BAMAObject<-function(from) { return(as(from, "BAMAObject")) }


# getStdCurv method
#   Computes the SC for a bamaobject
setGeneric("getStdCurv", function(object, ...) standardGeneric("getStdCurv"))
setMethod("getStdCurv", "BAMAObject", function(object, file)
{
  stdCtrls<-c(SB=0.0, S1=1.22, S2=4.88, S3=19.53, S4=78.13, S5=312.50, S6=1250.0, S7=5000.0) #pg/L
  tplate<-read.xls(file, skip=2, nrows=8)
  rownames(tplate)<-tplate[,1]
  tplate<-tplate[,2:13]
  #Cherche well_id - control
  #check exprs in ctrl
  ctrlCoords<-sapply(names(stdCtrls), function(x){
	coords<-which(tplate==x, arr.ind=TRUE)
	paste(LETTERS[coords[1]],coords[2],sep="")
	})
  ctrlFName<-sapply(ctrlCoords, function(x){
	pData(obj)[which(pData(obj)==x),"filename"]
	})

  return(res)
})
