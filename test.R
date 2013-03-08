library(LumiR)
library(ggplot2)
path <- "~/Dropbox/LENA P1/LENA_P1/"
system.time(blum<-read.experiment(path))
pData(blum)$control<-pData(blum)$well%in%paste0(LETTERS[1:7],1)
blum@phenoData@data$control<-pData(blum)$well%in%paste0(LETTERS[1:7],1)
blum@phenoData@data$concentration<-NA
blum@phenoData@data$concentration[pData(blum)$control]<-rev(c(1.22, 4.88, 19.53, 78.13, 312.50, 1250.0, 5000.0))
blum@phenoData@data$sample_type<-"human"
blum@phenoData@data$sample_type[pData(blum)$control]<-"standard"
blum@phenoData@data$sample_type[pData(blum)$well=="H1"]<-"blank"

df.bd<-melt(blum)
names(df.bd)[4]<-"Intensity"

df.bd.IL2<-subset(df.bd,analyte=="IL-2")
df.bd.IL2.standard<-subset(df.bd,analyte=="IL-2" & sample_type=="standard")

ggplot(df.bd.IL2,aes(x=well,y=Intensity))+geom_boxplot(aes(alpha=as.factor(concentration),fill=sample_type))+facet_wrap(~plate,nrow=2)+scale_y_log10()+theme_bw()+geom_boxplot(data=df.bd.IL2.standard,aes(alpha=as.factor(concentration)),fill="#4DAF4A")+scale_fill_brewer(palette="Set1")+ylab("Bead level intensity")+theme(text=element_text(size=16),axis.text.x=element_text(angle=55,hjust=1,size=8))+labs(fill="Sample type", alpha="Standard concentration")

ggsave("bead_level.png")

plot_layout(blum, fill="sample_type",alpha="as.factor(concentration)")+scale_fill_brewer(palette="Set1")+theme(text=element_text(size=16))+labs(fill="Sample type", alpha="Standard concentration")
ggsave("layout.png")


slum<-slummarize(blum)
df.mfi<-melt(slum)

df.mfi.subset<-subset(df.mfi,control==TRUE & analyte!="unknown0")
df.mfi.subset$analyte<-factor(df.mfi.subset$analyte)

ggplot(df.mfi.subset)+geom_point(aes(x=concentration,y=mfi,color=plate),alpha=.5,size=3)+scale_x_log10()+scale_y_log10()+facet_wrap(~analyte)+theme_bw()+geom_sc(object=slum,n=10,mapping=aes(x=concentration,y=mfi,color=plate),size=.8)+theme_bw()+scale_color_brewer(palette="Set1")+theme(text=element_text(size=16),panel.margin = unit(.1, "lines"))
ggsave("standard_curve.png")


+geom_boxplot(aes())



path<-"/Users/rgottard/Dropbox/Ofir's Luminex data/Listeria/Listeria/100816 wt and ActAFlaP time course/"

blum<-read.experiment(path)
bama.mfi<-BAMAsummarize(bama)
slum<-bSummarize(blum)
