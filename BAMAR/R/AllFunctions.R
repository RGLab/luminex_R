#----------
#MasterPlex

# read.Lxd
#   parse .lxd files to retrive assay/phenoData and MFI/count/mean/trimStuff
read.Lxd<-function(filename, path="./")
{
  lxd<-xmlTreeParse(filename)
  root<-xmlRoot(lxd)

  plate<-root[[which(names(root)=="Plate")]] #make a copy of the tree from the selected node
  platerows<-LETTERS[1:8]

  wVec<-c()
  lxbVec<-c()
  widx<-as.numeric(which(names(plate)=="Well"))
  wInfoVec<-vector("list", length(widx))
  for(curWidx in 1:length(widx))
  {
    #retrieve some well info
    ##TODO: parse <LocName> to put a name on each well
    curNode<-plate[[widx[[curWidx]]]]
    curWellInfo<-xmlAttrs(curNode)
    wellName<-paste(platerows[[as.integer(curWellInfo[["row"]])+1]], curWellInfo[["col"]], sep="")
    wVec<-c(wVec, wellName)
    curWellInfo<-c(Well=wellName, curWellInfo[3:length(curWellInfo)])
    wInfoVec[[curWidx]]<-curWellInfo
    #retrieve lxb filenames
    lxbStr<-xmlValue(curNode[["BinaryFileLoc"]])
    lxbFile<-strsplit(lxbStr, split="\\\\")[[1]]
    lxbFile<-lxbFile[length(lxbFile)]
    lxbVec<-c(lxbVec, lxbFile)
    #retrieve the beads
    bidx<-as.numeric(which(names(curNode)=="RSts"))
    #retrieve summary
    for(curBidx in 1:length(bidx))
    {
      curID<-as.numeric(xmlAttrs(curNode[[bidx[[curBidx]]]])[1])
      #curVals<-c(curWellInfo,xmlAttrs(curNode[[bidx[[curBidx]]]][4][[1]])) #4 is chan2 which should be the one used
      curSummary<-unlist(list(bid=curID,well_id=wellName,xmlAttrs(curNode[[bidx[[curBidx]]]][4][[1]])[c("median","mean","stdDev","cv")])) #4 is chan2 which should be the one used
      if(curWidx==1 & curBidx==1)
        summary<-data.frame(as.list(curSummary), stringsAsFactors=FALSE) #initialize
      else
        summary<-rbind(summary, as.list(curSummary)) 
    }
  }

  ##
  # phenoData
  ##
  pData<-data.frame(well_id=wVec, filename=lxbVec) #Minimal pD output: filename + well_id
  #pData<-cbind(t(data.frame(wInfoVec)), filename=lxbVec)
  #rownames(pData)<-seq(1:length(widx) #Maximal pD output: filename + all well info
  
  ##
  # assayData
  ## 
  setup<-root[[which(names(root)=="Setup")[1]]] #setup
  setupRun<-root[[which(names(root)=="Setup")[2]]] #setupRun
  ##TODO: which run is what?
  gateInfo<-xmlAttrs(setupRun[["SetGate"]])
  regionIdx<-as.numeric(which(names(setupRun)=="Region"))
  #initialize the data.frame with the control id
  aD<-data.frame(id=0, name="control", minCount=0, maxCount=100, active=1, stringsAsFactors=FALSE)
  #Parse other analytes
  for(curRIdx in 1:length(regionIdx))
  {
    curNode<-setupRun[[regionIdx[[curRIdx]]]]
    analyteInfo<-xmlAttrs(curNode)
    aD<-rbind(aD, as.list(analyteInfo))
  }
  
  ##
  # exprs
  ##
  #Create list with 1elt per well
  nFiles<-nrow(pData)
  exprsList<-vector('list', nFiles)
  cat("Looking for .lxb files in: ", path,"\n")
  cntFiles<-0
  for(fileIdx in 1:nFiles)
  {
    try(expr={
    lxbPath<-paste(path, pData[,"filename"][fileIdx], sep="/")
    #wID<-pData[,"well_id"][fileIdx]
    #I suppress the warnings about the int byte size as they are relevant only for the time variable that is not used
    suppressWarnings(lxb<-read.FCS(lxbPath))
    exprsLxb<-exprs(lxb)

    IDList<-sort(unique(exprsLxb[,"RID"]))
    sampleList<-vector('list', length(IDList)) #create list of set size
    for(idxID in 1:length(IDList))
    {   
      vals<-exprsLxb[which(exprsLxb[,"RID"]==IDList[idxID]), "RP1"]
      sampleList[[idxID]]<-vals
    }   
    names(sampleList)<-IDList
    exprsList[[fileIdx]]<-sampleList
    cntFiles<-cntFiles+1
    },#end expr
    silent=TRUE)
  }
  if(cntFiles==0)
  {
    warning("None of the binary files referenced in '",filename,"' were found in '",path,"'  The exprs slot will be empty.\nMake sure that 'path' leads to the folder where the .lxb files are located.\n")
  } else if(cntFiles<nFiles)
  {
    warning(nFiles-cntFiles, " files not found in '",path,"'\n")
  }
  names(exprsList)<-pData[,"well_id"]

   return(new("BAMAObject", phenoData=pData, assayData=aD, summary=summary, exprs=exprsList))
}


read.Lxb<-function(files)
{
  nFiles<-length(files)
  exprsList<-vector('list', nFiles)
  for(fileIdx in 1:nFiles)
  {
    #I suppress the warnings about the int byte size as they are relevant only for the time variable that is not used
    suppressWarnings(lxb<-read.FCS(files[[fileIdx]]))
    exprsLxb<-exprs(lxb)

    IDList<-sort(unique(exprsLxb[,"RID"]))
    sampleList<-vector('list', length(IDList)) #create list of set size
    for(idxID in 1:length(IDList))
    {   
      vals<-exprsLxb[which(exprsLxb[,"RID"]==IDList[idxID]), "RP1"]
      sampleList[[idxID]]<-vals
    }   
    names(sampleList)<-IDList
    exprsList[[fileIdx]]<-sampleList
    }
    names(exprsList)<-files
    return(exprsList)
}

#--------
# xPONENT
# Read raw data
.read.bd.xPONENT<-function(all.files,bid)
{
  #FIXME: There is a problem with the skip
  exprsList<-lapply(all.files,
                    function(x,bid){
                      con<-read.csv(x, skip=4, header=TRUE);
                      # Need to check the second argument and the number of lines to skip above
                      beadList<-split(con[,"RP1"],con[,2]);
                      bid.missing<-bid[bid%in%names(beadList)]
                      list.missing<-as.list(rep(NA,length(bid.missing)))
                      names(list.missing)<-bid.missing
                      beadList<-c(beadList,bid.missing)
                      # Order according to the bid
                      beadList<-beadList[order(names(beadList))]
                      return(beadList);},bid)
  
  return(exprsList)
}

# Function to combine raw and mapping data
read.luminex<-function(path="./")
{
  ## TODO add a function to read the file, sanitize it and check the format, column names, etc
  analyte.file<-list.files(path,pattern="analyte",full.names=TRUE)
  if(length(analyte.file)==0)
  {
    stop("Can't map bead ids to analyte, the 'analyte.csv' file is missing\n")
  }
  else
  {
    featureData<-as(read.csv(analyte.file,header=TRUE),'AnnotatedDataFrame')
  }

  ## TODO add a function to read the file, sanitize it and check the format, column names, etc
  ## Also check that the filenames match what's in the pheno file
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  if(length(pheno.file)==0)
  {
    warning("No pheno data provided, the 'phenotype.csv' file is missing\n")
  }
  else
  {    
  #To be added later, mapping file  
  #phenoData<-.read.mapping.pheno(pheno.file)
  }
  
  # Directory name (plate) only (without full path)
  plate.name<-list.dirs(path=path,recursive=F)
  plate.name<-unlist(lapply(plate.name,function(x)tail(strsplit(x,"/")[[1]],1)))
  # Make sure it's not a hidden directory
  which.to.keep<-sapply(plate.name,substr,1,1)!="."
  dirs<-list.dirs(path=path,recursive=F)[which.to.keep]
  plate.name<-plate.name[which.to.keep]
  
  # Get all files in all directories
  all.files<-unlist(lapply(dirs,list_files_with_exts,exts="csv"))
  filename<-unlist(lapply(all.files,function(x)tail(strsplit(x,"/")[[1]],1)))
  filename<-sub("^([^.]*).*", "\\1",filename)
  
  # Construct pheno data
  # Well id based on last 3 characters
  well<-unlist(lapply(filename,function(x){substr(x,nchar(x)-2,nchar(x))}))  
  # plate.name repeat unique plate name by number of wells
  plate.name<-rep(plate.name,sapply(plate.name,function(x,all.files){length(grep(x,all.files))},all.files))
  phenoData<-data.frame(plate=plate.name,filename=filename,well=well)
  phenoData<-as(phenoData,'AnnotatedDataFrame')
  
  # Read expression values
  exprs<-.read.bd.xPONENT(all.files,pData(featureData)$bid)
  
  names(exprs)<-pData(phenoData)$filename
  # Remove the outliers (bid not in the mapping file)        
  # Renames the analyte and replace the bid by the bead name
  exprs<-lapply(exprs,function(x,fd){x<-x[names(x)%in%fd$bid];names(x)<-fd$analyte[match(names(x),fd$bid)];return(x);},fd=pData(featureData))
  
  BAMAset<-new("BAMAset", phenoData=phenoData, featureData=featureData, exprs=exprs)
  
  return(BAMAset)
}

##TODO: Function to read our specified mapping format (csv or xls wld be OK)
#
#   INPUT: Read a mapping file
#   OUTPUT:
read.pheno.file<-function(file)
{
  pD<-read.csv(file)
}

### Summarize to MFIs and add standardCurves informations
BAMAsummarize<-function(from,type="MFI")
	  {
		  mat<-sapply(exprs(from),sapply,median)
		  mfiSet<-new("BAMAsummary", formula=as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f"))
		  exprs(mfiSet)<-mat
		  pData(mfiSet)<-pData(from)
		  fData(mfiSet)<-fData(from)
		  mfiSet@unit="MFI"

		  ##TODO: Potential args of the func
		  #hard coded 5-PL & inv
		  inv<-function(y, parmVec){exp(log(((parmVec[3] - parmVec[2])/(log(y) - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))}
		  
		  #fitInfo
		  df<-melt(mfiSet)
		  df<-subset(df, concentration!=0 & control==1)
		  nCtrl<-length(unique(df$well))
		  df.split<-split(df, df$analyte)
		  coeffs<-lapply(df.split, function(x){
					  res<-drm(log(mfi) ~ concentration, data=x,fct=LL.5())#, weights=1/mfi^2)
					  return(res$parmMat)
				  })
		  
		  calc_conc<-p100rec<-numeric(nrow(df))
		  for(idx in 1:nrow(df))
		  {
			  calc_conc[idx]<-inv(df[idx,"mfi"], coeffs[[df[idx,"analyte"]]])
			  p100rec[idx]<-calc_conc[idx]/df[idx,"concentration"]*100
		  }
		  li<-vector('list', 5)
		  names(li)<-c('b','c','d','e','f')
		  for(i in 1:5)
		  {
			  li[[i]]<-rep(sapply(coeffs, "[[", i), each=nCtrl)
		  }
		  df2<-cbind(df[,c("plate", "well", "analyte", "mfi", "concentration")], calc_conc, p100rec, li)
		  mfiSet@fit<-df2
		  
		mfiSet
	  }	  

	  
	  