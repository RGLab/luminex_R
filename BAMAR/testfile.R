
path <- "~/Dropbox/LENA P1/LENA_P1/"
bama<-read.luminex(path)
bama@phenoData@data$control<-pData(bama)$well%in%paste0(LETTERS[1:8],1)
bama@phenoData@data$concentration<-1
bama@phenoData@data$concentration[pData(bama)$control]<-rev(c(0.0, 1.22, 4.88, 19.53, 78.13, 312.50, 1250.0, 5000.0))

bama.mfi<-BAMAsummarize(bama)
df.mfi<-melt(bama.mfi)
df.mfi.subset<-subset(df.mfi,control==TRUE)

ggplot(df.mfi.subset)+geom_point(aes(x=concentration,y=mfi,color=plate),alpha=.5)+scale_x_log10()+scale_y_log10()+geom_sc(object=bama.mfi,n=10,mapping=aes(x=concentration,y=mfi,color=plate))+facet_wrap(~analyte)+theme_bw()

path<-"/Users/rgottard/Dropbox/Ofir's Luminex data/Listeria/Listeria/100816 wt and ActAFlaP time course"

bama<-read.luminex(path=path)
bama.mfi<-BAMAsummarize(bama)
