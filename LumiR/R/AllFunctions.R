#The root of the experiment
read.experiment<-function(path="./"){
  analyte.file<-list.files(path,pattern="analyte",full.names=TRUE)
  layout.file<-list.files(path,pattern="layout",full.names=TRUE)
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  plates<-list.dirs(path, recursive=FALSE)

  if(length(list_files_with_exts(plates[1], exts="lxb")>0)){
    type<-"LXB"; typeExt<-"lxb";
  } else if(length(list_files_with_exts(plates[1], exts="xml")>0)){
    type<-"BIOPLEX"; typeExt<-"xml";
  } else {
    type<-"XPONENT"; typeExt<-"csv"
  }
  
  all.files<-unlist(lapply(plates,list_files_with_exts,exts=typeExt))  

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
      len<-lapply(plates, function(x){length(list_files_with_exts(x,exts=typeExt))})
      fNames<-unlist(lapply(all.files, function(x)tail(strsplit(x,"/")[[1]],1)))
      plate<-rep(plates, len)
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

  if(length(layout.file)>0){
    layout<-.read.layout(layout.file)
    pData(phenoData)<-merge(pData(phenoData), layout, by="well")
  }

  #exprs
  if(type=="XPONENT"){exprs<-.read.exprs.xPonent(all.files)
  } else if(type=="LXB"){ exprs<-.read.exprs.lxb(all.files)
  } else {exprs<-.read.exprs.bioplex(all.files)}
  #names(exprs)<-plates

  #fData
  if(length(analyte.file)==0){#Only the bid available
    if(type=="LXB" & length(list_files_with_exts(path, exts="lxd"))>0){
      lxdFile<-list_files_with_exts(path, exts="lxd")
      featureData<-.read.lxd(lxdFile)
    } else{
      warning("Can't map bead ids to analyte, the 'analyte.csv' file is missing\n")
      bid<-unique(unlist(lapply(exprs, names)))
      analyte<-paste("unknown", bid, sep="")
      featureData<-as(data.frame(analyte=analyte, bid=bid), 'AnnotatedDataFrame')
    }
  } else { #checks: bid in mapping not in data and vice versa
    featureData<-.read.analyte(analyte.file)
  }
  beadInExprs<-unique(unlist(lapply(exprs, names)))
  notMappedBid<-beadInExprs[!beadInExprs%in%pData(featureData)$bid]
  if(length(notMappedBid>0)){
    cat("Removing",length(notMappedBid),"of the beads found in the data that are not found in the mapping file:", notMappedBid,"\n")
    exprs<-lapply(exprs, function(x){x[names(x)%in%pData(featureData)$bid]}) #remove non mapped beads (0=outliers)
    #pData(featureData)<-rbind(pData(featureData), data.frame(analyte=paste("unknown", notMappedBid, sep=""), bid=notMappedBid))
  }

  blum<-new("blum", phenoData=phenoData, featureData=featureData, exprs=exprs)
  return(blum)
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
    filenames<-unlist(lapply(filenames, function(x)paste(tail(strsplit(x,"/")[[1]],2), collapse="/")))
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
    slxb<-split(asdf$RP1,asdf$RID)
    exprsList[[fileIdx]]<-slxb
  }
  return(exprsList)
}
.read.exprs.bioplex<-function(filenames){
  exprs<-list()
  for(filename in filenames){
    xml<-xmlTreeParse(filename)
    root<-xmlRoot(xml)
    wNames<-xmlSApply(root[["Wells"]], function(x){
      paste(LETTERS[as.numeric(xmlAttrs(x)["RowNo"])], xmlAttrs(x)["ColNo"], sep="")
    })
    exprsFile<-xmlSApply(root[["Wells"]], function(x){
      str<-xmlValue(x[["BeadEventData"]])
      ss<-unlist(strsplit(str, split="\n"))
      sss<-strsplit(ss, split=" ")
      asdf<-as.data.frame(do.call(rbind, lapply(sss, as.numeric)))
      ret<-split(asdf[,2], asdf[,1])#bid
      return(ret)
    })
    names(exprsFile)<-paste(filename, wNames, sep="_")
    exprs<-c(exprs, exprsFile)
  }
  return(exprs)
}
  
  

.read.analyte<-function(analyte.file){
  df<-read.csv(analyte.file, header=TRUE)
  colnames(df)<-tolower(colnames(df))
  if(length(df)!=2 | !all(colnames(df)%in%c("analyte","bid"))) {
    stop("The analyte mapping file should be a csv file with two columns 'analyte' and 'bid'\n")
  } else {
    featureData<-as(df, "AnnotatedDataFrame")
  }
  return(featureData)
}

.read.layout<-function(layout.file){
  df<-read.csv(layout.file, header=TRUE)
  colnames(df)<-tolower(colnames(df))
  if(!all(colnames(df)%in%c("well", "sample_type", "concentration"))){#required cols
    stop("The layout mapping file should be a csv file with three columns 'well', 'sample_type' and 'concentration'\n")
  }
  if(length(unique(df$well))!=nrow(df)){#wells are unique
    stop("The layout file should contain only one line per well")
  }
  ##TODO: bkg=0 should be okay
  if(nrow(df[df$sample_type!="standard" & !is.na(df$concentration),])){#conc only set for standards
    stop("The 'concentration' in layout mapping file should only be set for standard wells\n Check wells: ",paste(as.character(df[df$sample_type!="standard" & !is.na(df$concentration),"well"]), collapse=","))
  }
  return(df)
}

.read.pheno.xPonent<-function(path, pheno.file){
  df <-read.csv(pheno.file, colClasses="factor")
  colnames(df)<-tolower(colnames(df))
  #df$concentration<-as.numeric(levels(df$concentration))[df$concentration] #More efficient than fact->char->numeric
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

### Summarize to MFIs and add standardCurves informations
slummarize<-function(from,type="MFI"){
  mat<-lapply(exprs(from),sapply,median)
  mat<-t(do.call("rbind",mat))
  mfiSet<-new("slum", formula=as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f"), inv=function(y, parmVec){exp(log(((parmVec[3] - parmVec[2])/(log(y) - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))}
  )
  exprs(mfiSet)<-mat
  pData(mfiSet)<-pData(from)
  fData(mfiSet)<-fData(from)
  mfiSet@unit="MFI"
  
  df<-melt(mfiSet)
  # subselects standards
  df<-subset(df, concentration!=0 & tolower(sample_type)=="standard")		  
  # Split by plate
  sdf<-split(df,df$plate)
  df2<-lapply(sdf,.fit_sc, mfiSet@inv)
  df2<-do.call("rbind",df2)
  mfiSet@fit<-df2

  pNames<-unlist(lapply(strsplit(colnames(mat), split="/"), "[[", 1))
  aNames<-fData(mfiSet)$analyte[match(rownames(mat), fData(mfiSet)$bid)]
  conc<-c()
  coefs<-unique(mfiSet@fit[,c("plate", "analyte", "b","c","d","e","f")])
  inv<-mfiSet@inv
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      conc<-c(conc, mfiSet@inv(mat[[i,j]], as.numeric(coefs[coefs$analyte==aNames[i] & coefs$plate==pNames[j],3:7])))
    }
  }
  concMat<-t(matrix(conc, ncol=ncol(mat)))
  assayData(mfiSet)<-list(exprs=mat, concentration=concMat)
  #mfiSet@assayData$concentration<-concMat
  mfiSet
}
	  

.fit_sc<-function(df, inv)
{
  
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
  sortCoeffs<-do.call("rbind", lapply(coeffs[df$analyte], t))
  colnames(sortCoeffs)<-c('b','c','d','e','f')
  df2<-cbind(df[,c("plate", "filename", "well", "analyte", "mfi", "concentration")], calc_conc, p100rec, sortCoeffs)
  return(df2)
}

results.conc.CSV<-function(object, file="./concentrations.csv"){
  mbs<-melt(object)
  concentration<-c()
  for(i in 1:nrow(mbs)){
    coefs<-getCoeffs(bs, mbs[i,"plate"], mbs[i, "analyte"])
    concentration<-c(concentration,as.numeric(inv(mbs[i, "mfi"], coefs)))
  }
  toWrite<-cbind(mbs[,c("plate", "well", "analyte", "mfi")], concentration)
  write.csv(toWrite, file=file, row.names=FALSE)
  
  return(invisible(toWrite))
}

results.curves.CSV<-function(object, file="./curves.csv"){
  bsfo<-bs@formula[3]
  fList<-c()
  bsfi<-unique(object@fit[,c("plate", "analyte", "b","c","d","e","f")])
  bsfi[,3:7]<-round(bsfi[,3:7], 4)
  for(i in 1:nrow(bsfi)){
    fList<-c(fList,gsub("c",bsfi[i,"c"],gsub("d",bsfi[i,"d"],gsub("e",bsfi[i,"e"],gsub("f",bsfi[i,"f"],bsfo)))))
  }
  toWrite<-cbind(sbs[,c("plate", "analyte")], Formula=fList)
  write.csv(toWrite, file=file, row.names=FALSE)
  return(invisible(toWrite))
}
