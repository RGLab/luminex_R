# Read the data and create a bamaset object
bama<-read.luminex(mapping="~/Dropbox/LENA P1/3-01-2012 LENA P1_20120301_145732.csv", path="~/Dropbox/LENA P1/3-01-2012 LENA P1_rcsv/")

# Add the concentration information
# Need to create a <- method for pData and fData
bama@phenoData@data$control<-pData(bama)$well%in%paste0(LETTERS[1:8],1)
bama@phenoData@data$concentration<-1
bama@phenoData@data$concentration[pData(bama)$control]<-rev(c(0.0, 1.22, 4.88, 19.53, 78.13, 312.50, 1250.0, 5000.0))

bsum<-BAMAsummarize(bama)
formula(bsum)<-getFormulas(bsum)
drawStdC(bsum)
p100<-getPercentRecov(bsum)

## Create dataframe for visualization
#df<-melt(bama)
#
## The object is big, so we subset first
#df2<-subset(df,analyte=="IL-2")
#
## Example visualization of raw data (slow!!)
#ggplot(df2)+geom_boxplot(aes(x=filename,y=intensity))
#
#
## Summarize (MFI)
#bama.mfi<-BAMAsummarize(bama)
#
## Create dataframe for visualization
#df.mfi<-melt(bama.mfi)
#
## Image type plot visualization
#ggplot(df.mfi)+geom_raster(aes(y=analyte,x=filename,fill=log(mfi)))+facet_wrap(~control)
#
#df.mfi.control<-subset(df.mfi,control==1)
#
#
### Create a dataframe of standard curves
#
#df.split<-split(df.mfi.control,df.mfi.control$analyte)
#df.sc<-lapply(df.split,function(x,n=100,weighted=FALSE,fct=LL.5())
#{
#  if(weighted)
#  {
#    weights=1/log(x$mfi)^2
#  }
#  else
#  {
#    weights<-NULL
#  }
#  concentration<-exp(seq(0,log(max(x$concentration)),length.out=n))
##  coef<-drm(mfi ~ concentration, data=x,fct=fct,weights=weights)$coefficients
##  b<-coef[1]
##  c<-coef[2]
##  d<-coef[3]
##  e<-coef[4]
##  f<-coef[5]  
##  
##  if(fct$name=="LL.5")
##    e<-log(e)
##  
##  mfi<-c + (d-c)/(1+exp(b*(log(concentration)-e)))^f
#  res<-drm(mfi ~ concentration, data=x,fct=fct,weights=weights)
#  mfi<-res$curve[[1]](concentration)
#  df<-data.frame(analyte=rep(x$analyte[1],n),concentration,mfi)
#  return(df)  
#},
#weighted=FALSE,fct=LL.5())
#df.sc<-as.data.frame(do.call("rbind",df.sc))
#
#ggplot(df.mfi.control)+geom_point(aes(y=mfi,x=concentration),size=2)+facet_wrap(~analyte)
#+scale_x_log10()+scale_y_log10()+theme_bw()+geom_line(data=df.sc,aes(y=mfi,x=concentration),color="blue")
#
#
##percentRecov
##Just use dfc mfi~conc and the inversion fct of the drm results
#res<-lapply(df.split,function(x)
#		{
#			res<-drm(mfi ~ concentration, data=x,fct=fct,weights=weights)
#		})
#
##For each row, do:
#
#
#
#recov<-numeric(nrow(df))
#calc_conc<-numeric(nrow(df))
#calc_conc<-numeric(15)
#recov<-numeric(15)
#for(i in 1:15)#nrow(dfc))
#{
#	calc_conc[i]<-res[[df[i,"analyte"]]]$fct$inversion(df[i,"mfi"], res[[df[i,"analyte"]]]$parmMat)
#	recov[i]<-calc_conc[i]/df[i,"concentration"]
#} ##Gotta remove the 0 first
#
#
#
#
#
#
#
#
#f5<-function(x,coef)
#{
#  b<-coef[1]
#  c<-coef[2]
#  d<-coef[3]
#  e<-coef[4]
#  f<-coef[5]  
#  return(c + (d-c)/(1+exp(b*(log(x)-log(e))))^f)
#}
#
#
## inverse = function (MFI,f, interval=c(0, 10000),...) {
#  solution<-uniroot(function (x,MFI,...){f(x,...)-MFI}, interval, MFI,...)[1]
#  return(solution)
#}
