# show method:
#   Shows class, slots, number of analytes, total number of measures
#
setMethod("show", "BAMAset", function(object){
  cat("An object of class BAMAset with",nrow(fData(object)),"analytes:","\n")
  cat("\t", as.character(head(fData(object)$analyte, 3)),"...", as.character(tail(fData(object)$analyte, 3)),"\n")
  cat(length(unlist(exprs(object), use.names=FALSE)), "measures of expression in", nrow(pData(object)),"wells.","\n")
  cat("And slots:", names(getSlots("BAMAset")),"\n")
})

#TODO: Add setters

# pData method:
#   phenoData accessor
setMethod("pData", "BAMAset", function(object){
  return(pData(object@phenoData))
})


## contains information about analytes
setMethod("fData", "BAMAset", function(object){
  return(pData(object@featureData))
})

# exprs accessor for bead level data a la eSet
setMethod("exprs", "BAMAset", function(object){
  return(object@exprs)
})

# fit accessor for standard curve fitting information
setGeneric("fit", function(object, ...) standardGeneric("fit"))
setMethod("fit", "BAMAsummary",function(object){
  return(object@fit)
})

# Subset method to subset a la eSet
setMethod("[","BAMAset",
          function(x,i,j,..., drop=FALSE)
          {
            if(!missing(i))
            {
              #Subset the samples
              bdata<-exprs(x)[j]
              #Subset the analytes
              bdata<-lapply(exprs(x),"[",i)              
            }
            else
            {
              #Subset the analytes         
              bdata<-lapply(exprs(x),"[",i)              
            }            
            newSet<-new('BAMAset'
                        ,exprs=bdata
                        ,phenoData=x@phenoData[j,]
                        ,featureData=x@featureData[i,])
            newSet            
          })

# I add the pheno and feature information by default we could add an option for this
# Need to sort this out with reshape2
# Also needs some work
setGeneric("melt",function(x,...){
  standardGeneric("melt")
})


setMethod("melt","BAMAset",
function(x)
  {
  # Use the melt function in reshape2
  # This generates a dataframe analyte, sample, intensity
  df<-reshape2::melt(exprs(x))
  names(df)<-c("intensity","analyte","filename")

#   l.sample<-length(exprs(x))
#   l.analyte<-length(exprs(x)[[1]])
#   l.bead.analyte<-sapply(exprs(x),function(x)sum(sapply(x,length)))
#   
#   pd.long<-sapply(1:l.sample,function(i,x,length)sapply(x[i,],rep,length[i]),length=l.bead.analyte,x=pData(x))                  
#   pd.long<-as.data.frame(do.call("rbind",pd.long))
#   
#   # Combine data and metadata
#   df<-cbind(df,pd.long)  

  
  df<-merge(df,pData(x),by="filename")
  df<-merge(df,fData(x),by="analyte")
  return(df)
  
})

setMethod("melt","BAMAsummary",
          function(x)
          {
            # Use the melt function in reshape2
            # This generates a dataframe analyte, sample, intensity
            df<-reshape2::melt(exprs(x))
            
            names(df)<-c("analyte","filename",tolower(x@unit))            
            #l.sample<-ncol(x)
            #l.analyte<-nrow(x)
            
            
#             pd.long<-apply(pData(x),1,rep,l.analyte)
#             pd.long<-lapply(1:l.sample,function(i,x,length)sapply(x[i,],rep,length),length=l.analyte,x=pData(x))
#             pd.long<-as.data.frame(do.call("rbind",pd.long))
#             
#             # Combine data and metadata
#             df<-cbind(df,pd.long)  
#            
            ## merge all information
            df<-merge(df,pData(x),by="filename")
            df<-merge(df,fData(x),by="analyte")
            return(df)
          })

#--------
# formula
# getter
#setGeneric("formula", function(object, ...) standardGeneric("formula"))
#setMethod("formula", "BAMAsummary", function(object){return(object@formula)})
#--------
# formula<-
# setter
setGeneric("formula<-", function(object, value, ...) standardGeneric("formula<-"))
setReplaceMethod("formula", "BAMAsummary", function(object, value){object@formula<-value; object})


setGeneric("geom_sc", function(object, n=100, data = NULL, 
				stat = "identity", position = "identity", 
				na.rm = FALSE, ...) standardGeneric("geom_sc"))

setMethod("geom_sc", "BAMAsummary",
          function(object, n=100, mapping = NULL,
                   stat = "identity", position = "identity", 
                   na.rm = FALSE, ...)
          {            
            
            fit<-object@fit
            # Extract the coefficients per plate/analyte
            # Split by plate then analyte                                
            sfit<-lapply(split(fit,fit$plate),function(x){split(x,x$analyte)})
            df.sc<-lapply(sfit,lapply,.df_sc_fit,object@formula)
            df.sc<-as.data.frame(do.call("rbind",lapply(df.sc,function(x)do.call("rbind",x))))
            
            ret<-geom_line(data=df.sc, mapping=mapping,
                           na.rm=na.rm, ...)
            return(ret)
          })

.df_sc_fit<-function(df,formula,n=100)
{
  # Concentration from min to max
  x<-exp(seq(0,log(max(df$concentration)),length.out=n))
  coef<-as.numeric(df[1,c('b','c','d','e','f')])
  # Should be able to parse this automatically
  b<-coef[1];c<-coef[2];d<-coef[3];e<-coef[4];f<-coef[5];
  # Evalulate formula
  mfi<-eval(parse(text=formula[3]))
  # Check whether the formula is fit on the log scale
  if(substr(as.character(formula[2]),1,3)=="log")
    mfi <- exp(mfi)  
  # basic dataframes with plate, filename, well, concentration, mfi
  newdf<-data.frame(plate=rep(df$plate[1],n),analyte=rep(df$analyte[1],n),concentration=x,mfi=mfi)
 return(newdf)
}
