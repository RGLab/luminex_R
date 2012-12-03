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
.read.bd.xPONENT<-function(path)
{
  files<-list_files_with_exts(path,"csv")

  exprsList<-lapply(files,function(x){con<-read.csv(x, skip=1, header=TRUE);beadList<-split(con[,"RP1"],con[,"RID"]);return(beadList);})
  return(exprsList)
}

# Read mapping
.read.mapping.xPONENT<-function(file)
{
  rl<-readLines(file, warn=FALSE)#A character vector 1char/line

  rl2<-gsub("\"", "", rl) #no quotes
  # Bead ID-Analyte mapping  
  position<-grep("DataType:,Units",rl2)
  analyte<-unlist(strsplit(rl2[position+1],","))[-1]
  bid<-unlist(strsplit(rl2[position+2],","))[-1]
  df<-data.frame(analyte=analyte,bid=bid)
  return(featureData=as(df,'AnnotatedDataFrame'))
}

# Function to combine raw and mapping data
read.luminex<-function(mapping.file=NULL, path="./")
{  
  if(is.null(mapping.file))
  {
    warning("You have not specified a mapping file, looking for data files in '",path,"'\n")
    ##TODO: try reading lxb in path and make a deefault aD/pD and calculate the summary
  }
  ext<-file_ext(mapping.file)
  if(ext=="lxd")
    return(read.Lxd(mapping.file, path))
  else if(ext=="csv")
  {
    featureData<-.read.mapping.xPONENT(mapping.file)
    # Read expression values
    exprs<-.read.bd.xPONENT(path)
    # Construct pheno data
    phenoData<-.makePhenoData.xPONENT(mapping.file,path)
    # Replace the names of the exprs object with the filename
    # If need we could construct a methods names for BAMAset to use better names
    
    names(exprs)<-pData(phenoData)$filename
    # Remove the outliers (bid not in the mapping file)        
    # Renames the analyte and replace the bid by the bead name
    exprs<-lapply(exprs,function(x,fd){x<-x[names(x)%in%fd$bid];names(x)<-fd$analyte[match(names(x),fd$bid)];return(x);},fd=pData(featureData))    
    BAMAset<-new("BAMAset", phenoData=phenoData, featureData=featureData, exprs=exprs)
  }
  return(BAMAset)
}

# .makePhenoData
#   INPUT: A character vector
#   OUTPUT: A data.frame that can be used for the pData slot
.makePhenoData.xPONENT<-function(mapping.file,path)
{
  files<-list_files_with_exts(path, "csv", all.files = FALSE, full.names=FALSE)
  # remove extension 
  files<-sub("^([^.]*).*", "\\1",files)
  # Plate information
  # Plate
  plate<-sub("^([^.]*).*", "\\1",mapping.file)
  plate<-tail(unlist(strsplit(plate,"/")),1)
  # Well (Assumes a standard file format)
  # Need to check it's always the case
  #well<-unlist(lapply(LETTERS[1:8],paste0,1:11))
  well<-unlist(lapply(files,function(x){tail(strsplit(x,"_")[[1]],1)}))
  
  phenoData<-data.frame(filename=files,plate=rep(plate,length(well)),well=well)
  return(as(phenoData,'AnnotatedDataFrame'))
}


##TODO: Function to read our specified mapping format (csv or xls wld be OK)
#
#   INPUT: Read a mapping file
#   OUTPUT:
read.pheno.file<-function(file)
{
  pD<-read.csv(file)
}


# getWellName
#   INPUT: int row, col;
#   OUTPUT: well_id
getWellName<-function(row, col)
{
  if(row < 1 | col < 1)
    stop("A well cannot have a NULL or negative row or column")
  if(row > 8 | col > 12)
    warning("The number of rows or column is too big for a standard 96 wells plate")
  return(paste(LETTERS[row], col, sep=""))
}

makeMapping<-function(folder)
{
  files<-list.files(folder)
  wells<-gsub(".csv", "", sapply(files,function(x){ul<-unlist(strsplit(x, split=c("_"))); return(ul[[length(ul)]])} ))
  cols<-as.numeric(substr(wells,2,nchar(wells)))
  rows<-substr(wells,1,1)
  rows<-as.numeric(sapply(rows, function(x) {which(LETTERS==x)} ))
  df<-data.frame(filename=files, row=rows, col=cols, row.names=wells)
  df <- df[with(df, order(row, col)), ]
  ##grp_idx: to match case with control
  ##type: case/ctrl/std/other ctrls
  ##dil: dilution for stdrds
  ##exp_conc: expected conc for stdrds
  ## Is it realistic to make defaults for that? Can I have a simple user input to specify that?:w
  
}


### Summarize to MFIs

BAMAsummarize<-function(from,type="MFI")
        {
        mat<-sapply(exprs(from),sapply,median)
        mfiSet<-new("BAMAsummary")
        exprs(mfiSet)<-mat
        pData(mfiSet)<-pData(from)
        fData(mfiSet)<-fData(from)
        mfiSet@unit="MFI"
        mfiSet
      }
