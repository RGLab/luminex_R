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

# Subset method to subset a la eSet
setMethod("[","BAMAset",
          function(x,i,j,..., drop=FALSE)
          {
            if(!missing(i))
            {
              #Subset the analytes
              bdata<-exprs(x)[i]
              #Subset the samples
              bdata<-lapply(exprs(x),"[",j)              
            }
            else
            {
              #Subset the samples              
              bdata<-lapply(exprs(x),"[",j)              
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
# getFormulas
# Returns the list of formulas for the BAMAsummary
setGeneric("getFormulas", function(object) standardGeneric("getFormulas"))
setMethod("getFormulas", "BAMAsummary",
		function(object)
		{
		  df<-melt(object)
		  browser()
		  df.ctrl<-subset(df, control==1)
		  df.split<-split(df.ctrl, df.ctrl$analyte) #split wwith analyte as factor
		  formulas.sc<-lapply(df.split, function(x){
					  res<-drm(mfi ~ concentration, data=x,fct=LL.5())
					  return(res$curve[[1]])
				  })
		  names(formulas.sc)<-names(df.split)
		  formula(object)<-formulas.sc
		})



#--------
# getPercentRecov
# Returns a d.f with mfi/conc/expected_conc/calculated_conc/percent_recovery BAMAsummary object
setGeneric("getPercentRecov", function(object, ...) standardGeneric("getPercentRecov"))
setMethod("getPercentRecov", "BAMAsummary",
		function(object)
		{
			fct<-LL.5()
#			weights=NULL
			if(length(formula(object))==0)
				stop("The formula slot is empty. Run 'getFormulas'to add the formulas for this object.\n")
			df<-melt(object)
			df<-subset(df, concentration!=0 & control==1, select=c("well", "analyte", "mfi", "concentration"))

			df.split<-split(df,df$analyte)
			res<-lapply(df.split,function(x)
					{
						weights<-1/log(x$mfi)^2
						res<-drm(mfi ~ concentration, data=x,fct=fct,weights=weights)
					})
			
			calc_conc<-numeric(nrow(df))
			recov<-numeric(nrow(df))
			for(idx in 1:nrow(df))
			{
				calc_conc[idx]<-res[[df[idx,"analyte"]]]$fct$inversion(df[idx,"mfi"], res[[df[idx,"analyte"]]]$parmMat)
				recov[idx]<-calc_conc[idx]/df[idx,"concentration"]
			}
			ret<-cbind(df, calc_conc, p100recov=recov)
			return(ret)
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


#--------
# drawStdC
# Returns the list of formulas for the BAMAsummary
setGeneric("ggplot_sc", function(object, ...) standadGeneric("drawStdC"))
setMethod("ggplot_sc", "BAMAsummary",
		function(object)
		{
			n<-100
			##TODO: get the formula from the formula...
			fct<-function(x, parmVec){parmVec[[2]] + (parmVec[[3]] - parmVec[[2]])/(1 + exp(parmVec[[1]] * (log(x) - log(parmVec[[4]]))))^parmVec[[5]]}

			fit<-object@fit
			concs<-exp(seq(0,log(max(fit$concentration)),length.out=n))
			coefs<-unique(fit[,c('b','c','d','e','f')])

			li<-lapply(unique(fit$analyte), function(x){
						fct(concs, coefs[x,])
					})
			stdf<-data.frame(analyte=rep(unique(fit$analyte), each=length(concs)), mfi=unlist(li), 
					concentration=rep(concs, length(unique(fit$analyte))))						
			
			ggplot(fit)+geom_point(aes(y=mfi,x=concentration),size=2)+facet_wrap(~analyte)+scale_x_log10()+scale_y_log10()+theme_bw()+geom_line(data=stdf,aes(y=mfi,x=concentration),color="blue")
			
		})


