
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
  well<-unlist(lapply(filename,function(x){gsub("_", "", substr(x,nchar(x)-2,nchar(x)))}))  
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
		  mat<-lapply(exprs(from),sapply,median)
      mat<-t(do.call("rbind",mat))
		  mfiSet<-new("BAMAsummary", formula=as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f"))
		  exprs(mfiSet)<-mat
		  pData(mfiSet)<-pData(from)
		  fData(mfiSet)<-fData(from)
		  mfiSet@unit="MFI"
		  
		  df<-melt(mfiSet)
      # subselects controls
		  df<-subset(df, concentration!=0 & control==1)		  
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
  
  nCtrl<-length(unique(df$well))
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
  return(df2)
}

	  
