#  #Bead level information
setClass("blum",
         representation=representation(
           ## Contains information about samples           
           phenoData="AnnotatedDataFrame",
           ## contains information about analytes           
           featureData="AnnotatedDataFrame",
           ## list of bead level data
           ## Stored as samples -> analytes
           exprs="list")
         )

setClass("bsum", 
		contains="ExpressionSet", 
		representation(unit="character",formula="formula", fit="data.frame"))


