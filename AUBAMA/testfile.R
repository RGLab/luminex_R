
#A is the MFI/RLU value for the minimum asymptote
#B is the Hill slope
#C is the concentration at the inflection point
#D is the MFI/RLU value for the maximum asymptote
#E is the asymmetry factor  #It should be F in drm

#Their params:
A<-11.53
B<--1.78
C<-1833.599
D<-8834.987
E<-0.813
para<-c(A=A,B=B,C=C,D=D,E=E)
paraNames<-c('A','B','C','D','E')
names(para)<-paraNames
func0<-function(x){para['A']+(para['D']/(1+(x/para['C'])^para['B'])^para['E'])}
plot(func0, from=1, log="xy", to=10000)

#data:
#in LENA_P1 folder
csv<-read.luminex(mapping="3-01-2012 LENA P1_20120301_145732.csv", path="3-01-2012 LENA P1_rcsv")

stdCtrls<-c(SB=0.0, S1=1.22, S2=4.88, S3=19.53, S4=78.13, S5=312.50, S6=1250.0, S7=5000.0) #pg/L
#Read the crappy assay_template in xls file
tplate<-read.xls("Assay_Template LENA P1 2-29-12.xls", skip=2, nrows=8)
rownames(tplate)<-tplate[,1]
tplate<-tplate[,2:13]
#Cherche well_id - control
ctrlCoords<-sapply(names(stdCtrls), function(x){
			coords<-which(tplate==x, arr.ind=TRUE)
			paste(LETTERS[coords[1]],coords[2],sep="")
		})
mfiLEPTIN<-MFI[2,ctrlCoords[2:8]]
dat<-data.frame(fi=mfiLEPTIN, expected_conc=stdCtrls[2:8])
datCalc<-c(5109.52, 1225.78, 333.25, 72.13, 21.82, 6.43, 1.21) #calculated conc with Holden's curve



#drc
resCalc1<-vector('list', 4) #from holden's calc conc to mfi
res<-vector('list', 4)
calc<-vector('list', 4) #from mfi to conc
#resCalc2<-vector('list', 4)
fcts<-list(l5(), L.5(), LL.5(), LL2.5())
names(resCalc1)<-names(resCalc2)<-c("l5", "L.5", "LL.5", "LL2.5")
#L.5 is without log(x)
#LL2.5 is with log(E)

for(i in 1:length(fcts))
{
  res[[i]]<-drm(fi ~ expected_conc, data=dat, fct=fcts[[i]])#, weights=dat$fi^-2)
  try({
    resCalc1[[i]]<-sapply(datCalc, function(x){res[[i]]$fct$bfct(x, res[[i]]$parmMat)})
  },silent=TRUE)
  calc[[i]]<-sapply(dat$fi, function(x){res[[i]]$fct$inversion(x, res[[i]]$parmMat)})
  plot(res[[i]], log="xy", main=names(resCalc1)[[i]])
}


#####################################################
#para<-res$coefficients
#names(para)<-paraNames
#func1<-function(x){para['C']+ ( (para['D']-para['C'])/
#				((1+exp((para['B']) * (log(x)-log(para['E']))))^para['F']) )}
#func2<-function(x){para['C']+ ((para['D']-para['C'])/
#				((1+exp((para['B']) * (log(x)-para['E'])))^para['F']) )} #This is when I am using LL2.5 since it already uses the log of E
#
#resCalc1[[i]]<-sapply(datCalc, func1)
#resCalc2[[i]]<-sapply(datCalc, func2)
