library(gdata)
library(multcomp)
library(drc)
library(Biobase)
library(ggplot2)
path<-"~/workspace/sexyTest/exp/"

bama<-read.luminex(path)
bama.mfi<-BAMAsummarize(bama)


df.mfi<-melt(sl)
df.mfi.subset<-subset(df.mfi,tolower(sample_type)=="standard")
ggplot(df.mfi.subset)+geom_point(aes(x=concentration,y=mfi,color=plate),alpha=.5)+scale_x_log10()+scale_y_log10()+geom_sc(object=sl)+facet_wrap(~analyte)+theme_bw()


sl<-slummarize(bl)
msl<-melt(sl)
ggplot2::ggplot(msl, ggplot2::aes(color=plate), alpha=0.5)+ggplot2::scale_x_log10()+ggplot2::scale_y_log10()+ggplot2::facet_wrap(~analyte)+geom_sc(sl)    +ggplot2::geom_point(ggplot2::aes(x=concentration, y=mfi))



k<-"TNF.A"
beads<-k
mSB<-sdSB<-list()
mSB[[k]]<-17.81107
sdSB[[k]]<-7.732063
deltas <- c(0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
PowerThresholds <- c(0.19,0.29,0.39,0.49,0.59,0.69,0.79)


#bama->mbama->split+mSB+sdSB->saxcyB->pVals+alpha

test<-function(){##TODO: This should prolly be saxcyB function that outputs the obj + pvals +alpha+deltas +thresh
	
	pList<-dList<-vector('list', length(PowerThresholds))
	names(pList)<-names(dList)<-as.character(PowerThresholds)
	for(k in beads){
		res2 <- saxcyBfit( san[[k]], mSB[[k]] )
		mfi2 <- getMFI2( res2$analyte )
		for ( fit in res2$fits ) {
			# select delta by power estimation
			power.est<-array(length(deltas))
			for ( delta.idx in 1:length(deltas) ) {
				power.est[delta.idx]<-mean(SAxCyBpower.ineq(fit, deltas[delta.idx], abs(coef(fit$ht)) ))
			}
			mfi <- drop.levels(subset(mfi2,control_idx==fit$ctrllvl) )
			an <- subset(res2$analyte,control_idx==fit$ctrllvl )
			# compact group ids
			an$group_name = drop.levels(an$group_name)
			#ctrl_name <- paste('control',fit$ctrllvl,sep='')
			ctrl_name <- as.character(unique(an[an$sample_type=="control", "group_name"]))
			#ctrl_name<-"wt_10^3_24h"
			# reset reference group id to control
			an$group_name<-relevel(drop.levels(an$group_name),ref=ctrl_name)
			
			qvals<-simplepval(fit, mfi, an )
			
			#save the delta and pvals
			for ( PowerThreshold in PowerThresholds ) {
				pt<-as.character(PowerThreshold)
				delta<-deltas[max(which(power.est>PowerThreshold))]
				dList[[pt]]<-c(dList[[pt]], rep(delta, length(levels(an$group_name))))
				pList[[pt]]<-c(pList[[pt]], SAxCyBpval(fit, delta))
			}
		}
	}
	return(list(pval=pList, delta=dList))
}#test

thresh<-selThreshold(dList)
pList[[as.character(thresh)]]
dList[[as.character(thresh)]]
