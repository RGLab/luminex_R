#The root of the experiment
read.experiment<-function(path="./"){
  analyte.file<-list.files(path,pattern="analyte",full.names=TRUE)
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  plates<-list.dirs(path, recursive=FALSE)

  if(length(list_files_with_exts(plates[1], exts="lxb")>0)){
    type<-"LXB"; typeExt<-"lxb";
  } else if(length(list_files_with_exts(plates[1], exts="xml")>0)){
    type<-"BIOPLEX"; typeExt<-"xml";
  } else {
    type<-"XPONENT"; typeExt<-"csv"
  }
  
  all.files<-lapply(plates,list_files_with_exts,exts=typeExt)
  
  #pData
  if(length(pheno.file)==0){#Get plate/file/well based on the structure of the folder
    warning("No pheno data provided, the 'phenotype.csv' file is missing\n")
    if(type=="BIOPLEX"){#The treatment should be different for BIOPLEX as there is only one file
      wells<-.getBioplexWellsID(all.files)
      wellsPerFile<-sapply(all.files, function(x){
        xmlSize(xmlRoot(xmlTreeParse(x))[["Wells"]])
      })
      fNames<-rep(unlist(lapply(all.files, lapply, function(x)tail(strsplit(x,"/")[[1]],1))), wellsPerFile)
      plate<-rep(plates, wellsPerFile)
      plate<-sapply(strsplit(plate, split="/"), tail, 1)
    } else {
      fNames<-unlist(lapply(all.files, lapply, function(x)tail(strsplit(x,"/")[[1]],1)))
      plate<-rep(plates, sapply(all.files, length))
      plate<-sapply(strsplit(plate, split="/"), tail, 1)
      if(type=="XPONENT"){
        wells<-.getXponentWellsID(fNames)
      } else if(type=="LXB"){
        wells<-.getLXBWellsID(fNames)
      } else {
        wells<-.getBioplexWellsID(fNames)
      }
    }
    phenoData<-data.frame(plate=plate,filename=fNames,well=wells)
    phenoData<-as(phenoData,'AnnotatedDataFrame')
  } else {
    phenoData<-.read.pheno.xPonent(path=path, pheno.file=pheno.file)
  }

  #exprs
  if(type=="XPONENT"){exprs<-lapply(all.files, .read.exprs.xPonent)
  } else if(type=="LXB"){ exprs<-lapply(all.files, .read.exprs.lxb)
  } else {exprs<-lapply(all.files, .read.exprs.bioplex)}
  names(exprs)<-plates

  #fData
  if(length(analyte.file)==0){#Only the bid available
    if(type=="LXB" & length(list_files_with_exts(path, exts="lxd"))>0){
      lxdFile<-list_files_with_exts(path, exts="lxd")
      featureData<-.read.lxd(lxdFile)
    } else{
      warning("Can't map bead ids to analyte, the 'analyte.csv' file is missing\n")
      bid<-unique(unlist(lapply(exprs, lapply, names)))
      analyte<-paste("unknown", bid, sep="")
      featureData<-as(data.frame(analyte=analyte, bid=bid), 'AnnotatedDataFrame')
    }
  } else { #checks: bid in mapping not in data and vice versa
    featureData<-.read.analyte.xPonent(analyte.file)
  }
  beadInExprs<-unique(unlist(lapply(exprs, lapply, names)))
  notMappedBid<-beadInExprs[!beadInExprs%in%pData(featureData)$bid]
  if(length(notMappedBid>0)){
    cat(length(notMappedBid),"of the beads found in the data are not found in the mapping file:", notMappedBid,"\n")
    pData(featureData)<-rbind(pData(featureData), data.frame(analyte=paste("unknown", notMappedBid, sep=""), bid=notMappedBid))
  }

  BAMAset<-new("BAMAset", phenoData=phenoData, featureData=featureData, exprs=exprs)
  return(BAMAset)
}

.getXponentWellsID<-function(filenames){
  filenames<-gsub(".csv", "", filenames)
  if(length(grep("Run", filenames))==0){#XPONENT v3.
    wellsID<-unlist(lapply(strsplit(filenames, split="_"), tail, 1))
  } else{#XPONENT v1.
    wellsNo<-as.numeric(gsub("Run", "", filenames))
    wellsID<-paste(LETTERS[(wellsNo-1)%%8+1],ceiling(wellsNo/8), sep="")
  }
  return(wellsID)
}
.getLXBWellsID<-function(filenames){#platename_wellID.lxb
  wellsID<-gsub(".lxb", "", unlist(lapply(strsplit(filenames, split="_"), tail, 1)))
  return(wellsID)
}
.getBioplexWellsID<-function(filenames){#file.xml > root>Wells
  wellsID<-unlist(lapply(filenames, function(x){
    root<-xmlRoot(xmlTreeParse(x))
    wNames<-xmlSApply(root[["Wells"]], function(x){
      paste(LETTERS[as.numeric(xmlAttrs(x)["RowNo"])], xmlAttrs(x)["ColNo"], sep="")
    })
    return(unlist(wNames))
  }))
  return(wellsID)
}

.read.exprs.xPonent<-function(filenames){
  sLine<-grep("[Ee]vent[Nn]o", readLines(filenames[1], n=5))-1
  exprsList<-lapply(filenames, function(x){
    con<-read.csv(x, skip=sLine, header=TRUE);
    beadList<-split(con[,"RP1"],con[,2]);
    beadList<-beadList[order(names(beadList))]
    return(beadList);})
    names(exprsList)<-filenames
  return(exprsList)
}
.read.exprs.lxb<-function(filenames){
  nFiles<-length(filenames)
  exprsList<-vector('list', nFiles)
  names(exprsList)<-filenames
  for(fileIdx in 1:nFiles){
    suppressWarnings(lxb<-read.FCS(filenames[[fileIdx]]))
    asdf<-as.data.frame(exprs(lxb))
    slxb<-split(asdf,asdf$RID)
    exprsList[[fileIdx]]<-slxb
  }
  return(exprsList)
}
.read.exprs.bioplex<-function(filenames){
  xml<-xmlTreeParse(filenames)
  root<-xmlRoot(xml)
  wNames<-xmlSApply(root[["Wells"]], function(x){
    paste(LETTERS[as.numeric(xmlAttrs(x)["RowNo"])], xmlAttrs(x)["ColNo"], sep="")
  })
  exprsList<-xmlSApply(root[["Wells"]], function(x){
    str<-xmlValue(x[["BeadEventData"]])
    ss<-unlist(strsplit(str, split="\n"))
    sss<-strsplit(ss, split=" ")
    asdf<-as.data.frame(do.call(rbind, lapply(sss, as.numeric)))
    ret<-split(asdf[,2], asdf[,1])#bid
    return(ret)
  })
  names(exprsList)<-wNames
  return(exprsList)
}
  

.read.analyte.xPonent<-function(analyte.file){
  df<-read.csv(analyte.file, header=TRUE)
  colnames(df)<-tolower(colnames(df))
  if(length(df)!=2 | !all(colnames(df)%in%c("analyte","bid"))) {
    stop("The analyte mapping file should be a csv file with two columns 'analyte' and 'bid'\n")
  } else {
    featureData<-as(df, "AnnotatedDataFrame")
  }
  return(featureData)
}

.read.pheno.xPonent<-function(path, pheno.file){
  df <-read.csv(pheno.file, colClasses="factor")
  colnames(df)<-tolower(colnames(df))
  df$concentration<-as.numeric(levels(df$concentration))[df$concentration] #More efficient than fact->char->numeric
  if(!all(c("plate","filename","well")%in%colnames(df))){
    stop("The phenotype mapping file must at least have the 'plate', 'filename' and 'well' columns\n")
  }
  for(i in 1:nrow(df)){
    if(length(list.files(paste(path, df[i,"plate"], sep="/"),pattern=as.character(df[i,"filename"])))==0){
      stop("The file ", as.character(df[i, "filename"]), " is not found in the given path. Verify plate and filename information in phenotype mapping file")
    }
  }
  phenoData<-as(df, "AnnotatedDataFrame")
  return(phenoData)
}

.read.lxd<-function(filename){#Workset>Setup>Region
  xml<-xmlTreeParse(filename)
  root<-xmlRoot(xml)
  setupNode<-root[["Setup"]]
  regionsIdx<-which(xmlSApply(setupNode, xmlName)=="Region")
  nRegion<-length(regionsIdx)
  ldf<-vector('list', nRegion)
  for(idx in 1:nRegion){
    ldf[[idx]]<-xmlAttrs(setupNode[[regionsIdx[[idx]]]])[c("name", "id")]
  }
  mat<-do.call(rbind, ldf)
  colnames(mat)<-c("analyte", "bid")
  featureData<-as(as.data.frame(mat), 'AnnotatedDataFrame')
  return(featureData)
}
