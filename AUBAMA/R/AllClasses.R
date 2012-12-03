#
#  Class to store information from BAMA experiment
#
setClass("BAMAset",
         representation=representation(
           ## Contains information about samples           
           phenoData="AnnotatedDataFrame",
           ## contains information about analytes           
           featureData="AnnotatedDataFrame",
           ## list of bead level data
           ## Stored as samples -> analytes
           exprs="list")
         )

setClass("BAMAsummary", contains="ExpressionSet", representation(unit="character",formula="list"))


