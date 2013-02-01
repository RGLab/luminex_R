
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
  sLine<-grep("[Ee]vent[Nn]o", readLines(all.files[1], n=5))-1
  exprsList<-lapply(all.files,
                    function(x,bid){
                      con<-read.csv(x, skip=sLine, header=TRUE);
                      # Need to check the second argument and the number of lines to skip above
                      beadList<-split(con[,"RP1"],con[,2]);
                      bid.missing<-bid[!bid%in%names(beadList)]
                      list.missing<-as.list(rep(NA,length(bid.missing)))
                      names(list.missing)<-bid.missing
                      beadList<-c(beadList,list.missing)
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
    #featureData<-as(read.csv(analyte.file,header=TRUE),'AnnotatedDataFrame')
    featureData<-.read.mapping.analyte(analyte.file)
  }


  ## TODO add a function to read the file, sanitize it and check the format, column names, etc
  ## Also check that the filenames match what's in the pheno file
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  if(length(pheno.file)==0)
  {
    warning("No pheno data provided, the 'phenotype.csv' file is missing\n")
    # Directory name (plate) only (without full path)
    plate.name<-list.dirs(path=path,recursive=FALSE)
    plate.name<-unlist(lapply(plate.name,function(x)tail(strsplit(x,"/")[[1]],1)))
    # Make sure it's not a hidden directory
    which.to.keep<-sapply(plate.name,substr,1,1)!="."
    dirs<-list.dirs(path=path,recursive=FALSE)[which.to.keep]
    plate.name<-plate.name[which.to.keep]
    # Get all files in all directories
    all.files<-unlist(lapply(dirs,list_files_with_exts,exts="csv"))
    filename<-unlist(lapply(all.files,function(x)tail(strsplit(x,"/")[[1]],1)))
    filename<-sub("^([^.]*).*", "\\1",filename)
    # Construct pheno data
    # Well id based on last 3 characters
    well<-unlist(lapply(filename,function(x){gsub("_", "", substr(x,nchar(x)-2,nchar(x)))}))  
    # plate.name repeat unique plate name by number of wells
    plate.name<-rep(plate.name,sapply(plate.name,function(x,all.files){length(grep(x,all.files))},all.files))
    phenoData<-data.frame(plate=plate.name,filename=filename,well=well)
    phenoData<-as(phenoData,'AnnotatedDataFrame')
  }
  else
  {    
    phenoData<-.read.mapping.pheno(pheno.file)
  }
  
  
  # Read expression values
  exprs<-.read.bd.xPONENT(paste(path,phenoData$plate,phenoData$filename,sep="/"),pData(featureData)$bid)
  
  names(exprs)<-pData(phenoData)$filename
  # Send warning regarding the outliers (bid not in the mapping file)        
  allBid<-unique(unlist(lapply(exprs, names)))
  notMappedBid<-allBid[!allBid%in%pData(featureData)$bid]
  if(length(notMappedBid>0))
  {
    cat(length(notMappedBid),"of the beads found in the data are not found in the mapping file:", notMappedBid,"\n")
  }
  pData(featureData)<-rbind(pData(featureData), data.frame(analyte=paste("unknown", notMappedBid, sep=""), bid=notMappedBid))
  # Renames the analyte and replace the bid by the bead name
  #exprs<-lapply(exprs,function(x,fd){x<-x[names(x)%in%fd$bid];names(x)<-fd$analyte[match(names(x),fd$bid)];return(x);},fd=pData(featureData))
  exprs<-lapply(exprs,function(x,fd){
    names(x)<-fd$analyte[match(names(x), fd$bid)]
    return(x)}, fd=pData(featureData))
  #exprs<-lapply(exprs,function(x,fd){
    #newNames<-as.character(fd$analyte[match(names(x), fd$bid)])
    #newNames[is.na(newNames)]<-names(x)[is.na(newNames)]
    #names(x)<-newNames
    #return(x)}, fd=pData(featureData))
  
  BAMAset<-new("BAMAset", phenoData=phenoData, featureData=featureData, exprs=exprs)
  
  return(BAMAset)
}

#INPUT: filename
#OUTPUT AnnotatedDataFrame object
.read.mapping.analyte<-function(analyte.file)
{
  df<-read.csv(analyte.file, header=TRUE)
  colnames(df)<-tolower(colnames(df))
  if(length(df)!=2 | !all(colnames(df)%in%c("analyte","bid")))
  {
    stop("The analyte mapping file should be a csv file with two columns 'analyte' and 'bid'\n")
  }
  else
  {
    featureData<-as(df, "AnnotatedDataFrame")
  }
  return(featureData)
}

#INPUT: filename
#OUTPUT: AnnotatedDataFrame object
.read.mapping.pheno<-function(pheno.file)
{
  df <-read.csv(pheno.file, colClasses="factor")
  colnames(df)<-tolower(colnames(df))
  phenoData<-as(df, "AnnotatedDataFrame")
  return(phenoData)
}

### Summarize to MFIs and add standardCurves informations
BAMAsummarize<-function(from,type="MFI")
	  {
		  mat<-lapply(exprs(from),sapply,median)
      mat<-t(do.call("rbind",mat))
		  mfiSet<-new("BAMAsummary", formula=as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f"))
		  exprs(mfiSet)<-mat
		  pData(mfiSet)<-pData(from)
		  fData(mfiSet)<-fData(from)
		  mfiSet@unit="MFI"
		  
		  df<-melt(mfiSet)
      # subselects standards
		  df<-subset(df, concentration!=0 & tolower(sample_type)=="standard")		  
		  # Split by plate
		  sdf<-split(df,df$plate)
      df2<-lapply(sdf,.fit_sc)
      df2<-do.call("rbind",df2)
      
		  mfiSet@fit<-df2
		  
		mfiSet
	  }	  

.fit_sc<-function(df)
{
  inv<-function(y, parmVec){exp(log(((parmVec[3] - parmVec[2])/(log(y) - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))}
  
  nCtrl<-length(unique(df$well)) #number of wells with standards
  df.split<-split(df, df$analyte)
  coeffs<-lapply(df.split, function(x){
    res<-drm(log(mfi) ~ concentration, data=x,fct=LL.5())
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
  df2<-cbind(df[,c("plate", "filename", "well", "analyte", "mfi", "concentration")], calc_conc, p100rec, li)
  print(nCtrl)
  return(df2)
}

read.raw.bioplex<-function(filename){
  xml<-xmlTreeParse(filename)
  root<-xmlRoot(xml)
  
}

.read.well<-function(node){
  infos<-xmlAttrs(node) #Events: number of measures. NumberColumns. Column_1, Column_2.. : colnames.
  str<-xmlValue(node) #a single str with " " and  "\n"
  ss<-unlist(strsplit(str, split="\n"))
  sss<-strsplit(ss, split=" ")
  df<-do.call(rbind, lapply(sss, as.numeric))
  colnames(df)<-infos[-c(1,2)] ##TODO: probably useless, should be done after merging everything
  return(df)
}  
