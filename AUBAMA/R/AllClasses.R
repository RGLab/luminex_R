#
#  Class to store information from BAMA experiment
#
setClass("BAMAObject",
        representation=representation(
                phenoData="data.frame",
                assayData="data.frame",
		summary="data.frame",
                exprs="list") #a list of list mapped into pData
        )
