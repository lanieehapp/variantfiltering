##single sample VCF merge functions


library(vcfR)
library(stringr)
library(plyr)
library(parallel)
library(doParallel)
library(data.table)
library(bedr)



three.caller.merge<-function(samp, file.list){
  
  all.data<-read.vcfR(file=file.list[1])
  all.hc.data<-parse.hc.vcf(all.data)
  
  all.data<-read.vcfR(file=file.list[2])
  all.s2.data<-parse.s2.vcf(all.data)
  
  all.data<-read.vcfR(file=file.list[3])
  all.dv.data<-parse.dv.vcf(all.data)
  
  all.data.merged<-merge.vcfs(samp,all.hc.data, all.dv.data, all.s2.data)
  
  
  
  gt.summary<-get.gt.summary(all.data.merged)
  all.fix.merged<-data.frame(all.data.merged[[1]], stringsAsFactors = FALSE)
  all.fix.merged<-all.fix.merged[order(all.fix.merged$CHROM_POS_REF_ALT),]
  all.info.merged<-all.data.merged[[2]]
  all.gt.merged<-all.data.merged[[3]]
  
  print("After initial merge:")
  print("fix.filt:")
  print(head(all.fix.merged$CHROM_POS_REF_ALT))
  print("info.filt:")
  print(head(all.info.merged$CHROM_POS_REF_ALT))
  print("gt.filt:")
  print(head(all.gt.merged$CHROM_POS_REF_ALT))
  
  all.gt.merged<-data.frame(all.gt.merged[,1],gt.summary, all.gt.merged[,2:ncol(all.gt.merged)], stringsAsFactors = FALSE, check.names = FALSE)
  colnames(all.gt.merged)[1]<-"CHROM_POS_REF_ALT"
  
  ##remove deep variant double entry variants
  if(nrow(all.fix.merged) != nrow(all.gt.merged)){
    tmp<-data.frame(table(all.gt.merged$CHROM_POS_REF_ALT), stringsAsFactors = FALSE)
    bad.vars<-as.character(tmp[tmp$Freq>1,1])
  }
  if(exists("bad.vars")){
    for( var in bad.vars){
      all.info.merged<-all.info.merged[!grepl(var, all.info.merged$CHROM_POS_REF_ALT),]
      all.gt.merged<-all.gt.merged[!grepl(var, all.gt.merged$CHROM_POS_REF_ALT),]
      all.fix.merged<-all.fix.merged[!grepl(var, all.fix.merged$CHROM_POS_REF_ALT),]
    }
  }
  
  all.info.merged<-all.info.merged[!grepl("ERCC|GL|KI", all.info.merged$CHROM_POS_REF_ALT),]
  all.gt.merged<-all.gt.merged[!grepl("ERCC|GL|KI", all.gt.merged$CHROM_POS_REF_ALT),]
  all.fix.merged<-all.fix.merged[!grepl("ERCC|GL|KI", all.fix.merged$CHROM),]
  
  all.fix.merged$REF<-as.character(all.fix.merged$REF)
  all.fix.merged$ALT<-as.character(all.fix.merged$ALT)
  all.info.merged<-cbind(all.info.merged, get.annovar.filters(all.fix.merged))
  
  ##apply basic filter and caller filter
  all.fix.filt<-all.fix.merged[all.info.merged$Caller_Filter & all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
  all.gt.filt<-all.gt.merged[all.info.merged$Caller_Filter & all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
  all.info.filt<-all.info.merged[all.info.merged$Caller_Filter & all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
  
  print("After basic and caller filters:")
  print("fix.filt:")
  print(head(all.fix.filt$CHROM_POS_REF_ALT))
  print("info.filt:")
  print(head(all.info.filt$CHROM_POS_REF_ALT))
  print("gt.filt:")
  print(head(all.gt.filt$CHROM_POS_REF_ALT))
  
  ##calculate filter for good, rare variants in windows 
  all.vars<-NULL
  
  for(chr in unique(all.fix.filt$CHROM)){
    
    vars<-all.fix.filt[all.fix.filt$CHROM==chr,]
    vars$POS<-as.numeric(vars$POS)
    vars<-vars[order(vars$POS),]
    
    vars$CLUST<-NA
    if(nrow(vars)==1){
      vars$CLUST<-FALSE
    }
    if(nrow(vars)>1){
      vars$CLUST[1]<-(as.numeric(vars$POS[2])-as.numeric(vars$POS[1])<10)
      for(i in 2:(nrow(vars))){
        vars$CLUST[i]<-((as.numeric(vars$POS[i])-as.numeric(vars$POS[i-1])<10) | (as.numeric(vars$POS[i+1])-as.numeric(vars$POS[i])<10))
      }
      end<-nrow(vars)
      vars$CLUST[end]<-(as.numeric(vars$POS[end])-as.numeric(vars$POS[end-1])<10)
    }
    
    
    all.vars<-rbind(all.vars, vars)  
  }
  
  print(all.vars[1:5, 1:5])
  
  all.vars<-all.vars[order(all.vars$CHROM_POS_REF_ALT),]
  all.info.filt$CLUST<-all.vars$CLUST
  #colnames(all.info.merged)[ncol(all.info.merged)]<-paste0(samp, ".CLUST")
  
  ##AF filter
  all.info.filt$AF_Filter<-all.gt.filt$afMax<0.25
  
  ##depth filter
  all.info.filt$Depth_Filter<-all.gt.filt$dpMax>9
  
  
  all.fix.filt<-all.fix.filt[!(all.info.filt$CLUST & all.info.filt$AF_Filter) & all.info.filt$Depth_Filter ,]
  all.gt.filt<-all.gt.filt[!(all.info.filt$CLUST & all.info.filt$AF_Filter) & all.info.filt$Depth_Filter ,]
  all.info.filt<-all.info.filt[!(all.info.filt$CLUST & all.info.filt$AF_Filter) & all.info.filt$Depth_Filter ,]
  
  # Remove any rows with NA as the position index CHROM_POS_REF_ALT
  print(sprintf("Rows with NA in CHROM_POS_REF_ALT in all.fix.filt: %d/%d",
    sum(is.na(all.fix.filt$CHROM_POS_REF_ALT)), nrow(all.fix.filt)))
  print(sprintf("Rows with NA in CHROM_POS_REF_ALT in all.info.filt: %d/%d", 
    sum(is.na(all.info.filt$CHROM_POS_REF_ALT)), nrow(all.info.filt)))
  print(sprintf("Rows with NA in CHROM_POS_REF_ALT in all.gt.filt: %d/%d", 
    sum(is.na(all.gt.filt$CHROM_POS_REF_ALT)), nrow(all.gt.filt)))
  
  rows_to_include <- which(!is.na(all.fix.filt$CHROM_POS_REF_ALT))
  all.fix.filt = all.fix.filt[rows_to_include,]
  all.info.filt = all.info.filt[rows_to_include,]
  all.gt.filt = all.gt.filt[rows_to_include,]
  
  print("Sizes after removing NAs:")
  print(sprintf("Rows with NA in CHROM_POS_REF_ALT in all.fix.filt: %d/%d",
    sum(is.na(all.fix.filt$CHROM_POS_REF_ALT)), nrow(all.fix.filt)))
  print(sprintf("Rows with NA in CHROM_POS_REF_ALT in all.info.filt: %d/%d", 
    sum(is.na(all.info.filt$CHROM_POS_REF_ALT)), nrow(all.info.filt)))
  print(sprintf("Rows with NA in CHROM_POS_REF_ALT in all.gt.filt: %d/%d", 
    sum(is.na(all.gt.filt$CHROM_POS_REF_ALT)), nrow(all.gt.filt)))
    
  all.fix.wl<-all.fix.merged[all.info.merged$whitelist,]
  all.info.wl<-all.info.merged[all.info.merged$whitelist,]
  all.gt.wl<-all.gt.merged[all.info.merged$whitelist,]
  
  print("Before returning data, after applying depth and cluster filter:")
  print("fix.filt:")
  print(head(all.fix.filt$CHROM_POS_REF_ALT))
  print("info.filt:")
  print(head(all.info.filt$CHROM_POS_REF_ALT))
  print("gt.filt:")
  print(head(all.gt.filt$CHROM_POS_REF_ALT))
  
  #colnames(all.gt.merged)[2:ncol(all.gt.merged)]<-paste(samp, colnames(all.gt.merged)[2:ncol(all.gt.merged)], sep=".")
  #colnames(all.info.merged)[2:ncol(all.info.merged)]<-paste(samp, colnames(all.info.merged)[2:ncol(all.info.merged)], sep=".")
  
  #colnames(all.gt.filt)[2:ncol(all.gt.filt)]<-paste(samp, colnames(all.gt.filt)[2:ncol(all.gt.filt)], sep=".")
  #colnames(all.info.filt)[2:ncol(all.info.filt)]<-paste(samp, colnames(all.info.filt)[2:ncol(all.info.filt)], sep=".")
  
  #colnames(all.gt.wl)[2:ncol(all.gt.wl)]<-paste(samp, colnames(all.gt.wl)[2:ncol(all.gt.wl)], sep=".")
  #colnames(all.info.wl)[2:ncol(all.info.wl)]<-paste(samp, colnames(all.info.wl)[2:ncol(all.info.wl)], sep=".")
  
  
  return(list(all.fix.merged, all.info.merged, all.gt.merged, all.fix.filt, all.info.filt, all.gt.filt, all.fix.wl, all.info.wl,all.gt.wl))
  
}

parse.hc.vcf<-function(all.data){
  all.fix<-data.frame(all.data@fix, stringsAsFactors = FALSE)
  CHROM_POS_REF_ALT<-paste(all.fix$CHROM, all.fix$POS, all.fix$REF, all.fix$ALT, sep="-")
  all.gt<-data.frame(all.data@gt, stringsAsFactors = FALSE)
  
  all.info<-all.fix$INFO
  
  print("Parsing INFO column")
  tmp<-mclapply(all.info, function(x){
    
    num.values<-str_count(x, ";")
    parsed<-str_split_fixed(x, ";", n=num.values+1)
    
    colnames(parsed)<-unlist(lapply(parsed, function(x){strsplit(x, "=")[[1]][1]}))
    parsed<-gsub(".*=", "", parsed)
    parsed<-data.frame(parsed, stringsAsFactors = FALSE)
    return(parsed)
  })
  
  all.info.parsed<-rbind.fill(tmp)
  if(!("ANNOVAR_DATE" %in% colnames(all.info.parsed))){
    all.info.hc<-cbind(all.info.parsed, all.fix[,6:7])
    colnames(all.info.hc)<-paste0("HC_", colnames(all.info.hc))
    all.fix.hc<-all.fix[,1:5]
  }
  
  if("ANNOVAR_DATE" %in% colnames(all.info.parsed)){
    annovar.start<-which(colnames(all.info.parsed)=="ANNOVAR_DATE")
    annovar.end<-which(colnames(all.info.parsed)=="ALLELE_END")
    if(annovar.end == ncol(all.info.parsed)){
      all.info.annovar<-all.info.parsed[,annovar.start:ncol(all.info.parsed)]
      all.info.hc<-cbind(all.info.parsed[,1:(annovar.start-1)], all.fix[,6:7])
    }
    if(annovar.end != ncol(all.info.parsed)){
      all.info.annovar<-all.info.parsed[,annovar.start:annovar.end]
      all.info.hc<-cbind(all.info.parsed[,1:(annovar.start-1)], all.info.parsed[,(annovar.end+1):ncol(all.info.parsed)])
      all.info.hc<-cbind(all.info.hc, all.fix[,6:7])
    }
    
    colnames(all.info.hc)<-paste0("HC_", colnames(all.info.hc))
    all.fix.hc<-cbind(all.fix[,1:5], all.info.annovar)
  }
  
  all.fix.hc<-all.fix.hc[,!(colnames(all.fix.hc) %in% c("CIGAR", "RU", "REFREP", "IDREP"))]
  
  print(paste("Parsed ", ncol(all.info.parsed), " INFO columns."))
  
  
  print("Parsing GTs")
  tmp<-mclapply(1:nrow(all.gt), function(x){
    sel<-all.gt[x,]
    num.values<-str_count(sel[2], ":")
    parsed<-as.vector(t(str_split_fixed(sel[1,2:ncol(sel)], ":", n=num.values+1)))
    col.type<-str_split_fixed(sel[1], ":", n=num.values+1)  
    parsed<-as.data.frame(t(parsed), stringsAsFactors = FALSE)
    new.names<-paste("HC", col.type, sep="_")
    colnames(parsed)<-new.names
    
    parsed<-data.frame(parsed, stringsAsFactors = FALSE)
    
    return(parsed)
  })
  
  all.gt.parsed.hc<-rbind.fill(tmp)
  print(paste("Parsed ", ncol(all.gt.parsed.hc), " GT columns."))
  
  nt<-c("A", "T", "C", "G")
  for(i in 1:nrow(all.info.hc)){
    if (all.fix.hc$REF[i] %in% nt & all.fix.hc$ALT[i] %in% nt){
      all.info.hc$HC_FILTER[i]<-(as.numeric(all.info.hc[i,"HC_QD"])>2 & !(is.na(all.info.hc[i, "HC_QD"])))& (as.numeric(all.info.hc[i,"HC_MQ"])>40& !(is.na(all.info.hc[i, "HC_MQ"])))& (as.numeric(all.info.hc[i,"HC_FS"])<60 & !(is.na(all.info.hc[i, "HC_FS"]))) & (as.numeric(all.info.hc[i,"HC_SOR"])<3 & !(is.na(all.info.hc[i, "HC_SOR"]))) & (as.numeric(all.info.hc[i,"HC_MQRankSum"])>(-12.5)& !(is.na(all.info.hc[i, "HC_MQRankSum"])))& (as.numeric(all.info.hc[i,"HC_ReadPosRankSum"])>(-8)& !(is.na(all.info.hc[i, "HC_ReadPosRankSum"])))
    }
    if (!(all.fix.hc$REF[i] %in% nt & all.fix.hc$ALT[i] %in% nt)){
      all.info.hc$HC_FILTER[i]<-(as.numeric(all.info.hc[i,"HC_QD"])>2 & !(is.na(all.info.hc[i, "HC_QD"]))) & (as.numeric(all.info.hc[i,"HC_FS"])<200 & !(is.na(all.info.hc[i, "HC_FS"]))) & (as.numeric(all.info.hc[i,"HC_SOR"])<10 & !(is.na(all.info.hc[i, "HC_SOR"]))) & (as.numeric(all.info.hc[i,"HC_ReadPosRankSum"])>(-20)& !(is.na(all.info.hc[i, "HC_ReadPosRankSum"])))
    }
  }
  
  #all.fix.hc$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  #all.info.hc$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  #all.gt.parsed.hc$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  
  all.fix.hc<-cbind(CHROM_POS_REF_ALT, all.fix.hc)
  all.info.hc<-cbind(CHROM_POS_REF_ALT, all.info.hc)
  all.gt.parsed.hc<-cbind(CHROM_POS_REF_ALT, all.gt.parsed.hc)
  
  
  return(list(all.fix.hc, all.info.hc, all.gt.parsed.hc)) 
}

parse.s2.vcf<-function(all.data){
  all.fix<-data.frame(all.data@fix, stringsAsFactors = FALSE)
  CHROM_POS_REF_ALT<-paste(all.fix$CHROM, all.fix$POS, all.fix$REF, all.fix$ALT, sep="-")
  all.gt<-data.frame(all.data@gt, stringsAsFactors = FALSE)
  
  all.info<-all.fix$INFO
  
  print("Parsing INFO column")
  tmp<-mclapply(all.info, function(x){
    
    num.values<-str_count(x, ";")
    parsed<-str_split_fixed(x, ";", n=num.values+1)
    
    colnames(parsed)<-unlist(lapply(parsed, function(x){strsplit(x, "=")[[1]][1]}))
    parsed<-gsub(".*=", "", parsed)
    parsed<-data.frame(parsed, stringsAsFactors = FALSE)
    return(parsed)
  })
  
  all.info.parsed<-rbind.fill(tmp)
  if(!("ANNOVAR_DATE" %in% colnames(all.info.parsed))){
    all.info.s2<-cbind(all.info.parsed, all.fix[,6:7])
    colnames(all.info.s2)<-paste0("S2_", colnames(all.info.s2))
    all.fix.s2<-all.fix[,1:5]
  }
  
  if("ANNOVAR_DATE" %in% colnames(all.info.parsed)){
    annovar.index<-which(colnames(all.info.parsed)=="ANNOVAR_DATE")
    all.info.annovar<-all.info.parsed[,annovar.index:ncol(all.info.parsed)]
    all.info.s2<-cbind(all.info.parsed[,1:(annovar.index-1)], all.fix[,6:7])
    colnames(all.info.s2)<-paste0("S2_", colnames(all.info.s2))
    all.fix.s2<-cbind(all.fix[,1:5], all.info.annovar)
  }
  
  
  all.fix.s2<-all.fix.s2[,!(colnames(all.fix.s2) %in% c("CIGAR", "RU", "REFREP", "IDREP"))]
  
  print(paste("Parsed ", ncol(all.info.parsed), " INFO columns."))
  
  
  print("Parsing GTs")
  tmp<-mclapply(1:nrow(all.gt), function(x){
    sel<-all.gt[x,]
    num.values<-str_count(sel[2], ":")
    parsed<-as.vector(t(str_split_fixed(sel[1,2:ncol(sel)], ":", n=num.values+1)))
    col.type<-str_split_fixed(sel[1], ":", n=num.values+1)  
    parsed<-as.data.frame(t(parsed), stringsAsFactors = FALSE)
    new.names<-paste( "S2", col.type, sep="_")
    colnames(parsed)<-new.names
    
    parsed<-data.frame(parsed, stringsAsFactors = FALSE)
    
    return(parsed)
  })
  
  all.gt.parsed.s2<-rbind.fill(tmp)
  print(paste("Parsed ", ncol(all.gt.parsed.s2), " GT columns."))
  
  
  
  all.fix.s2$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  all.info.s2$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  all.gt.parsed.s2$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  
  return(list(all.fix.s2, all.info.s2, all.gt.parsed.s2))
}

parse.dv.vcf<-function(all.data){
  all.fix<-data.frame(all.data@fix, stringsAsFactors = FALSE)
  CHROM_POS_REF_ALT<-paste(all.fix$CHROM, all.fix$POS, all.fix$REF, all.fix$ALT, sep="-")
  all.gt<-data.frame(all.data@gt, stringsAsFactors = FALSE)
  
  all.info<-all.fix$INFO
  
  print("Parsing INFO column")
  tmp<-mclapply(all.info, function(x){
    
    num.values<-str_count(x, ";")
    parsed<-str_split_fixed(x, ";", n=num.values+1)
    
    colnames(parsed)<-unlist(lapply(parsed, function(x){strsplit(x, "=")[[1]][1]}))
    parsed<-gsub(".*=", "", parsed)
    parsed<-data.frame(parsed, stringsAsFactors = FALSE)
    return(parsed)
  })
  
  all.info.parsed<-rbind.fill(tmp)
  
  if(!("ANNOVAR_DATE" %in% colnames(all.info.parsed))){
    all.info.dv<-cbind(all.info.parsed, all.fix[,6:7])
    colnames(all.info.dv)<-paste0("DV_", colnames(all.info.dv))
    all.fix.dv<-all.fix[,1:5]
  }
  if("ANNOVAR_DATE" %in% colnames(all.info.parsed)){
    annovar.index<-which(colnames(all.info.parsed)=="ANNOVAR_DATE")
    all.info.annovar<-all.info.parsed[,annovar.index:ncol(all.info.parsed)]
    all.info.parsed<-all.info.parsed[,1:annovar.index-1]
    all.info.dv<-cbind(all.info.parsed, all.fix[,6:7])
    colnames(all.info.dv)<-paste0("DV_", colnames(all.info.dv))
    all.fix.dv<-cbind(all.fix[,1:5], all.info.annovar)
  }
  
  
  all.fix.dv<-all.fix.dv[,!(colnames(all.fix.dv) %in% c("CIGAR", "RU", "REFREP", "IDREP"))]
  
  
  print("Parsing GTs")
  tmp<-mclapply(1:nrow(all.gt), function(x){
    sel<-all.gt[x,]
    num.values<-str_count(sel[2], ":")
    parsed<-as.vector(t(str_split_fixed(sel[1,2:ncol(sel)], ":", n=num.values+1)))
    col.type<-str_split_fixed(sel[1], ":", n=num.values+1)  
    parsed<-as.data.frame(t(parsed), stringsAsFactors = FALSE)
    new.names<-paste("DV", col.type, sep="_")
    colnames(parsed)<-new.names
    
    parsed<-data.frame(parsed, stringsAsFactors = FALSE)
    
    return(parsed)
    
  })
  
  all.gt.parsed.dv<-rbind.fill(tmp)
  print(paste("Parsed ", ncol(all.gt.parsed.dv), " GT columns."))
  
  all.fix.dv$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  all.info.dv$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  all.gt.parsed.dv$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  
  
  return(list(all.fix.dv, all.info.dv, all.gt.parsed.dv))
}

merge.vcfs<-function(samp, all.hc.data, all.dv.data, all.s2.data){
  
  #merged.fix<-unique(rbind(all.hc.data[[1]], all.dv.data[[1]], all.s2.data[[1]]))
  merged.fix<-unique(rbind.fill(all.hc.data[[1]], all.dv.data[[1]], all.s2.data[[1]]))
  #merged.fix = merged.fix[order(merged.fix$CHROM_POS_REF_ALT),]
  merged.info<-Reduce(function(x,y) merge(x=x, y=y, by="CHROM_POS_REF_ALT", all.x=T, all.y=T), list(all.hc.data[[2]], all.dv.data[[2]], all.s2.data[[2]]))
  
  merged.gt<-Reduce(function(x,y) merge(x=x, y=y, by="CHROM_POS_REF_ALT", all.x=T, all.y=T), list(all.hc.data[[3]], all.dv.data[[3]], all.s2.data[[3]]))
  
  
  merged.info$Caller_Filter<-(merged.info$`HC_FILTER`=="TRUE" & !(is.na(merged.info$`HC_FILTER`))) | (merged.info$`DV_FILTER`=="PASS"& !(is.na(merged.info$`DV_FILTER`))) | (merged.info$`S2_FILTER`=="PASS"& !(is.na(merged.info$`S2_FILTER`)))
  
  merged.fix<-merged.fix[!duplicated(merged.fix$CHROM_POS_REF_ALT),]
  
  
  return(list(merged.fix, merged.info, merged.gt))
  
}

get.gt.summary<-function(all.data.merged){
  
  
  all.gt.parsed<-all.data.merged[[3]]
  CHROM_POS_REF_ALT<-all.gt.parsed$CHROM_POS_REF_ALT
  samp<-colnames(all.gt.parsed)[2]
  samp<-unlist(strsplit(samp, "[.]"))[1]
  sel.ad<-all.gt.parsed[,grep(".AD$", colnames(all.gt.parsed))]
  
  ncallers<-rowSums(!(is.na(sel.ad)))
  
  af.calc<-apply(sel.ad, 1, function(x)str_split(x, ","))
  
  #outer outer list for each variant
  #outer list for each caller
  #inner list of REF[1] and ALT[2] reads
  
  tmp<-suppressWarnings(lapply(af.calc, function(x){
    lapply(x, function(y){
      if(is.na(y)){
        return(0)
      }
      if(!(is.na(y))){
        return(as.numeric(y[2])/(as.numeric(y[1])+as.numeric(y[2])))
      }
    })
  }))
  af.max<-lapply(tmp, function(x){
    return(max(unlist(x)))
  })
  
  sel.dp<-as.data.frame(sapply(all.gt.parsed[,grep(".DP", colnames(all.gt.parsed))], as.numeric))
  sel.dp[is.na(sel.dp)]<-0
  
  dp.max<-apply(sel.dp, 1, max)
  
  
  
  samp.meta<-data.frame(ncallers,unlist(af.max), dp.max, stringsAsFactors = FALSE)
  
  colnames(samp.meta)<-c("nCallers", "afMax", "dpMax")
  return(samp.meta)
  
}

get.annovar.filters<-function(all.fix.merged){
  
  pop.freq.data<-suppressWarnings(sapply(all.fix.merged[,grep("gnomad|ExAC|esp6500|X1000g",colnames(all.fix.merged),value=T)], as.numeric))
  pop.freq.data[is.na(pop.freq.data)] <- 0
  all.fix.merged$pop.freq.max.all <- apply(pop.freq.data,1,max)
  rm(pop.freq.data)
  
  all.filters<-data.frame(all.fix.merged$CHROM_POS_REF_ALT, stringsAsFactors = FALSE)
  all.filters$exonic <- all.fix.merged$Func.refGene=="exonic"
  load("filtered_whitelist_09302020.RData")
  #load("variantFilter/filtered_whitelist_09302020.RData")
  all.filters$whitelist<-all.fix.merged$CHROM_POS_REF_ALT %in% wl$CHROM_POS_REF_ALT
  all.filters$rare.variants <-all.fix.merged$pop.freq.max.all<=0.01
  all.filters$not.syn <- all.fix.merged[,'ExonicFunc.refGene']!="synonymous_SNV" 
  all.filters$not.superdups <- all.fix.merged$genomicSuperDups=="."
  all.filters$basic.filters<- (all.filters$whitelist) | (all.filters$rare & all.filters$not.superdups)  
  all.filters$pop.freq.max.all<-all.fix.merged$pop.freq.max.all
  
  
  ###repeatMasker
  
  load("repeat_masker.RData")
  options(scipen = 999)
  
  print(head(all.fix.merged$REF))
  
  end.pos<-unlist(lapply(1:nrow(all.fix.merged), function(x){as.numeric(paste(all.fix.merged$POS[x])) + max(nchar(as.character(all.fix.merged$REF[x])), nchar(as.character(all.fix.merged$ALT[x])))}))
  
  var.locs<-c(paste0(all.fix.merged$CHROM, ":", all.fix.merged$POS,"-", end.pos))
  #var.locs<-cbind(all.fix.merged[,1:2], all.fix.merged[,2]+1)
  table(is.valid.region(var.locs))
  var.locs<-bedr.sort.region(var.locs)
  
  
  not.repeatmasker<-!(in.region(var.locs, bed.file.filt))
  
  names(not.repeatmasker)<-var.locs
  not.repeatmasker<-not.repeatmasker[order(names(not.repeatmasker))]
  all.filters$not.repeatmasker<-not.repeatmasker
  
  
  return(all.filters[,2:9])
  
}



get_whitelist_vars_dna<-function(dna_bam, ref){
  whitelist_path<-"filtered_whitelist_09302020.txt"
  
  comm<-paste0('parallel --colsep "\t" samtools mpileup -a -l ', whitelist_path, ' --fasta-ref ',ref, ' ', dna_bam, ' -r {1} :::: ', ref_fai ,' > dna_whitelist.txt' )
  
  system(comm)
  
  dna_whitelist<-read.csv(file="dna_whitelist.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE,  quote="")
  colnames(dna_whitelist)<-c("Chr", "Pos", "Ref", "Depth", "mpileup", "Qual")
  
  dna_whitelist$ref_FS<-str_count(dna_whitelist$mpileup, pattern="[.]")
  dna_whitelist$ref_RS<-str_count(dna_whitelist$mpileup, pattern=",")
  dna_whitelist$A_FS<-str_count(dna_whitelist$mpileup, pattern="A")
  dna_whitelist$A_RS<-str_count(dna_whitelist$mpileup, pattern="a")
  dna_whitelist$C_FS<-str_count(dna_whitelist$mpileup, pattern="C")
  dna_whitelist$C_RS<-str_count(dna_whitelist$mpileup, pattern="c")
  dna_whitelist$G_FS<-str_count(dna_whitelist$mpileup, pattern="G")
  dna_whitelist$G_RS<-str_count(dna_whitelist$mpileup, pattern="g")
  dna_whitelist$T_FS<-str_count(dna_whitelist$mpileup, pattern="T")
  dna_whitelist$T_RS<-str_count(dna_whitelist$mpileup, pattern="t")
  dna_whitelist$DEL<-str_count(dna_whitelist$mpileup, pattern="[*]")
  dna_whitelist$INS<-str_count(dna_whitelist$mpileup, pattern="[+]")
  
  dna_whitelist$A_evidence<-dna_whitelist$A_FS>0 & dna_whitelist$A_RS>0
  dna_whitelist$C_evidence<-dna_whitelist$C_FS>0 & dna_whitelist$C_RS>0
  dna_whitelist$T_evidence<-dna_whitelist$T_FS>0 & dna_whitelist$T_RS>0
  dna_whitelist$G_evidence<-dna_whitelist$G_FS>0 & dna_whitelist$G_RS>0
  
  tmp<-dna_whitelist[dna_whitelist$A_evidence,]
  tmp$Alt_depth<-tmp$A_FS+tmp$A_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  
  
  vars<-cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("A", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF)
  
  tmp<-dna_whitelist[dna_whitelist$C_evidence,]
  tmp$Alt_depth<-tmp$C_FS+tmp$C_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  vars<-rbind(vars, cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("C", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF))
  
  tmp<-dna_whitelist[dna_whitelist$T_evidence,]
  tmp$Alt_depth<-tmp$T_FS+tmp$T_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  vars<-rbind(vars, cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("T", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF))
  
  tmp<-dna_whitelist[dna_whitelist$G_evidence,]
  tmp$Alt_depth<-tmp$G_FS+tmp$G_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  vars<-data.frame(rbind(vars, cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("G", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF)), stringsAsFactors = FALSE)
  
  
  colnames(vars)<-c("Chr", "Pos", "Ref", "Alt", "DNA_depth_total", "DNA_depth_alt", "DNA_AF")
  vars$CHROM_POS_REF_ALT<-paste(vars$Chr, vars$Pos, vars$Ref, vars$Alt, sep="-")
  
  return(vars)
  
}

get_whitelist_vars_rna<-function(rna_bam, ref){
  whitelist_path<-"filtered_whitelist_09302020.txt"
  comm<-paste0('parallel --colsep "\t" samtools mpileup -a -l ', whitelist_path, ' --fasta-ref ',ref, ' ', rna_bam, ' -r {1} :::: ', ref_fai ,' > rna_whitelist.txt'  )
  
  system(comm)
  
  rna_whitelist<-read.csv(file="rna_whitelist.txt", sep="\t", stringsAsFactors = FALSE, header = FALSE, quote="")
  
  colnames(rna_whitelist)<-c("Chr", "Pos", "Ref", "Depth", "mpileup", "Qual")
  
  rna_whitelist$ref_FS<-str_count(rna_whitelist$mpileup, pattern="[.]")
  rna_whitelist$ref_RS<-str_count(rna_whitelist$mpileup, pattern=",")
  rna_whitelist$A_FS<-str_count(rna_whitelist$mpileup, pattern="A")
  rna_whitelist$A_RS<-str_count(rna_whitelist$mpileup, pattern="a")
  rna_whitelist$C_FS<-str_count(rna_whitelist$mpileup, pattern="C")
  rna_whitelist$C_RS<-str_count(rna_whitelist$mpileup, pattern="c")
  rna_whitelist$G_FS<-str_count(rna_whitelist$mpileup, pattern="G")
  rna_whitelist$G_RS<-str_count(rna_whitelist$mpileup, pattern="g")
  rna_whitelist$T_FS<-str_count(rna_whitelist$mpileup, pattern="T")
  rna_whitelist$T_RS<-str_count(rna_whitelist$mpileup, pattern="t")
  rna_whitelist$DEL<-str_count(rna_whitelist$mpileup, pattern="[*]")
  rna_whitelist$INS<-str_count(rna_whitelist$mpileup, pattern="[+]")
  
  rna_whitelist$A_evidence<-rna_whitelist$A_FS>0 & rna_whitelist$A_RS>0
  rna_whitelist$C_evidence<-rna_whitelist$C_FS>0 & rna_whitelist$C_RS>0
  rna_whitelist$T_evidence<-rna_whitelist$T_FS>0 & rna_whitelist$T_RS>0
  rna_whitelist$G_evidence<-rna_whitelist$G_FS>0 & rna_whitelist$G_RS>0
  
  tmp<-rna_whitelist[rna_whitelist$A_evidence,]
  tmp$Alt_depth<-tmp$A_FS+tmp$A_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  
  
  rna_vars<-cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("A", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF)
  
  tmp<-rna_whitelist[rna_whitelist$C_evidence,]
  tmp$Alt_depth<-tmp$C_FS+tmp$C_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  rna_vars<-rbind(rna_vars, cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("C", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF))
  
  tmp<-rna_whitelist[rna_whitelist$T_evidence,]
  tmp$Alt_depth<-tmp$T_FS+tmp$T_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  rna_vars<-rbind(rna_vars, cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("T", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF))
  
  tmp<-rna_whitelist[rna_whitelist$G_evidence,]
  tmp$Alt_depth<-tmp$G_FS+tmp$G_RS
  tmp$AF<-tmp$Alt_depth/tmp$Depth
  rna_vars<-data.frame(rbind(rna_vars, cbind(tmp$Chr, tmp$Pos, tmp$Ref, rep("G", nrow(tmp)), tmp$Depth, tmp$Alt_depth, tmp$AF)), stringsAsFactors = FALSE)
  
  colnames(rna_vars)<-c("Chr", "Pos", "Ref", "Alt", "RNA_depth_total", "RNA_depth_alt", "RNA_AF")
  
  rna_vars$CHROM_POS_REF_ALT<-paste(rna_vars$Chr, rna_vars$Pos, rna_vars$Ref, rna_vars$Alt, sep="-")
  
  return(rna_vars)
}

merge_whitelists<-function(dna_bam, rna_bam, ref, single.sample.merged){
  
  dna_whitelist<-get_whitelist_vars_dna(dna_bam, ref)
  rna_whitelist<-get_whitelist_vars_rna(rna_bam, ref)
  
  merged_whitelist<-merge(dna_whitelist, rna_whitelist, by="CHROM_POS_REF_ALT", all.x=TRUE, all.y=TRUE)
  merged_whitelist<-cbind(merged_whitelist$CHROM_POS_REF_ALT, merged_whitelist[,grep("DNA|RNA", colnames(merged_whitelist))])
  colnames(merged_whitelist)[1]<-"CHROM_POS_REF_ALT"
  #colnames(merged_whitelist)[2:ncol(merged_whitelist)]<-paste0(samp, ".",colnames(merged_whitelist)[2:ncol(merged_whitelist)])
  
  load("filtered_whitelist_09302020.RData")
  
  all_whitelist_annot<-merge(wl, merged_whitelist, by="CHROM_POS_REF_ALT", all.x=FALSE, all.y=FALSE)
  
  all_whitelist_discovery<-cbind(single.sample.merged[[7]], single.sample.merged[[8]], single.sample.merged[[9]])
  CHROM_POS_REF_ALT<-all_whitelist_discovery$CHROM_POS_REF_ALT
  all_whitelist_discovery<-cbind(CHROM_POS_REF_ALT, all_whitelist_discovery[,!(grepl("CHROM_POS", colnames(all_whitelist_discovery)))])
  
  all_whitelist_merged<-merge(all_whitelist_discovery, all_whitelist_annot, by="CHROM_POS_REF_ALT", all.x=TRUE, all.y=TRUE)
  
  #handle case where discovery had no whitelist variants
  if (nrow(all_whitelist_merged)>0 & ncol(all_whitelist_merged[,(grepl("[.]x", colnames(all_whitelist_merged)))])>0 & ncol(all_whitelist_merged[,(grepl("[.]y", colnames(all_whitelist_merged)))])>0){
    all_whitelist_x<-all_whitelist_merged[,(grepl("[.]x", colnames(all_whitelist_merged)))]
    all_whitelist_y<-all_whitelist_merged[,(grepl("[.]y", colnames(all_whitelist_merged)))]
    all_whitelist_other<-all_whitelist_merged[,!(grepl("[.]x|[.]y", colnames(all_whitelist_merged)))]
    
    for(i in 1:nrow(all_whitelist_merged)){
      if(is.na(all_whitelist_x[i,1])){
        all_whitelist_x[i,]<-all_whitelist_y[i,]
      }
    }
    
    all_whitelist_merged<-cbind(all_whitelist_x, all_whitelist_other)
  }
  
  #handle case where no whitelist variants exist
  if(nrow(all_whitelist_merged)==0){
    print("No whitelist variants found")
    load("whitelist_column_names.RData")
    
    all_whitelist_merged<-data.frame(matrix(NA, nrow=0, ncol=207))
    colnames(all_whitelist_merged)<-wl_cols
    return(all_whitelist_merged)
  }
  
  
  all_whitelist_a<-all_whitelist_merged[,grepl("CHROM_POS_REF_ALT", colnames(all_whitelist_merged))]
  all_whitelist_b<-all_whitelist_merged[,!(grepl("CHROM_POS_REF_ALT", colnames(all_whitelist_merged)))]
  
  all_whitelist_merged<-cbind(all_whitelist_a, all_whitelist_b)
  
  #consolidated chrom,pos,ref,alt from merging
  bad_cols<-c("CHROM", "POS", "REF", "ALT", "Chr", "Start", "End", "Ref", "Alt")
  
  all_whitelist_merged<-all_whitelist_merged[,!(colnames(all_whitelist_merged) %in% bad_cols)]
  
  colnames(all_whitelist_merged)[1]<-"CHROM_POS_REF_ALT"
  
  tmp<-unlist(strsplit(all_whitelist_merged[,1], "-"))
  tmp<-t(data.frame(lapply(all_whitelist_merged[,1], function(x) {unlist(strsplit(x, "-"))})))
  
  all_whitelist_merged<-cbind(tmp, all_whitelist_merged)
  
  colnames(all_whitelist_merged)[1:4]<-c("CHROM", "POS", "REF", "ALT")
  
  return(all_whitelist_merged) 
}

get_validation_vars_rna<-function(rna_bam, ref, single.sample.merged, samp){
  
  filt.fix<-single.sample.merged[[4]]
  
  print("input filt.fix:")
  print(head(filt.fix))

  end_pos<-unlist(lapply(1:nrow(filt.fix), function(x){as.numeric(filt.fix$POS[x]) + max(nchar(filt.fix$REF[x]), nchar(filt.fix$ALT[x]))}))  

  filt.bed<-cbind(filt.fix$CHROM, as.numeric(filt.fix$POS)-1, end_pos)

  filt.bed<-unique(filt.bed)

  filt.bed.file <- paste0(getwd(), "/filt_bed.tmp.txt")

  write.table(filt.bed, sep="\t", file=filt.bed.file, row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  comm<-paste0('parallel --colsep "\t" samtools mpileup -a -l ', filt.bed.file, ' --fasta-ref ',ref, ' ', rna_bam, ' -r {1} :::: ', ref_fai ,' > rna_filt.txt' )
  
  system(comm)
  
  rna_filt<-read.csv(file="rna_filt.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE, quote = "", fill=FALSE)
  
  colnames(rna_filt)<-c("Chr", "Pos", "Ref", "Depth", "mpileup", "Qual")
  rna_filt$CUR_POS<-paste(rna_filt$Chr, rna_filt$Pos, sep="-")
  
  filt.fix$CHROM_POS_REF_ALT<-paste(filt.fix$CHROM, filt.fix$POS, filt.fix$REF, filt.fix$ALT, sep="-")
  
  rna_filt$ref_FS<-str_count(rna_filt$mpileup, pattern="[.]")
  rna_filt$ref_RS<-str_count(rna_filt$mpileup, pattern=",")
  rna_filt$A_FS<-str_count(rna_filt$mpileup, pattern="A")
  rna_filt$A_RS<-str_count(rna_filt$mpileup, pattern="a")
  rna_filt$C_FS<-str_count(rna_filt$mpileup, pattern="C")
  rna_filt$C_RS<-str_count(rna_filt$mpileup, pattern="c")
  rna_filt$G_FS<-str_count(rna_filt$mpileup, pattern="G")
  rna_filt$G_RS<-str_count(rna_filt$mpileup, pattern="g")
  rna_filt$T_FS<-str_count(rna_filt$mpileup, pattern="T")
  rna_filt$T_RS<-str_count(rna_filt$mpileup, pattern="t")
  rna_filt$DEL<-str_count(rna_filt$mpileup, pattern="[*]")
  rna_filt$INS<-str_count(rna_filt$mpileup, pattern="[+]")
  
  rna_filt$A_evidence<-rna_filt$A_FS>0 & rna_filt$A_RS>0
  rna_filt$C_evidence<-rna_filt$C_FS>0 & rna_filt$C_RS>0
  rna_filt$T_evidence<-rna_filt$T_FS>0 & rna_filt$T_RS>0
  rna_filt$G_evidence<-rna_filt$G_FS>0 & rna_filt$G_RS>0
  rna_filt$INS_evidence<-rna_filt$INS>0
  rna_filt$DEL_evidence<-rna_filt$DEL>0
  
  nt<-c("A", "C", "T", "G")
  all_rna<-NULL
  for(i in 1:nrow(filt.fix)){
    curr.pos<-paste(filt.fix$CHROM[i], filt.fix$POS[i], sep="-")
    curr.ref<-filt.fix$REF[i]
    curr.alt<-filt.fix$ALT[i]
    tmp<-rna_filt[rna_filt$CUR_POS == curr.pos,]
    
    RNA_depth_total<-tmp$Depth
    if(RNA_depth_total>=30){
      RNA_evidence="Covered"
    } else if(RNA_depth_total < 2){
      RNA_evidence="No coverage"
    } else if(RNA_depth_total<30){
      RNA_evidence<-"Low coverage"
    }
    
    #processing for SNVs
    if(curr.alt %in% nt & curr.ref %in% nt){
      if(sum(tmp$A_evidence, tmp$T_evidence, tmp$C_evidence, tmp$G_evidence, tmp$INS_evidence, tmp$DEL_evidence)>1){
        RNA_evidence<-paste(RNA_evidence, "Multiallelic Locus", sep=";")
      }
      
      if(curr.alt == "A" & tmp$A_evidence){
        RNA_depth_alt<-tmp$A_FS+tmp$A_RS
        RNA_AF<-RNA_depth_alt/RNA_depth_total
        RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
      } else if(curr.alt == "C" & tmp$C_evidence){
        RNA_depth_alt<-tmp$C_FS+tmp$C_RS
        RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
        RNA_AF<-RNA_depth_alt/RNA_depth_total
      } else if(curr.alt == "G" & tmp$G_evidence){
        RNA_depth_alt<-tmp$G_FS+tmp$G_RS
        RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
        RNA_AF<-RNA_depth_alt/RNA_depth_total
      } else if(curr.alt == "T" & tmp$T_evidence){
        RNA_depth_alt<-tmp$T_FS+tmp$T_RS
        RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
        RNA_AF<-RNA_depth_alt/RNA_depth_total
      } else {
        RNA_depth_alt=0
        RNA_evidence<-paste(RNA_evidence, "FALSE", sep=";")
        RNA_AF<-0
      }
      
    }
    
    #processing for indels
    if(!(curr.alt %in% nt) | !(curr.ref %in% nt) ){
      RNA_depth_alt=0
      start_pos<-filt.fix$POS[i]
      end_pos<-max(nchar(curr.alt), nchar(curr.ref))+as.numeric(filt.fix$POS[i])-1
      
      all_pos<-unlist(lapply(start_pos:end_pos, function(x){paste(filt.fix$CHROM[i], x, sep="-")}))
      for(pos in all_pos){
        tmp2<-rna_filt[rna_filt$CUR_POS == pos,]
        if(tmp2$DEL_evidence){
          RNA_depth_alt=max(RNA_depth_alt, tmp2$DEL)
          RNA_depth_total<-max(RNA_depth_total, tmp2$Depth)
          RNA_evidence<-paste(RNA_evidence, "DEL", sep=";")
        } else if (tmp2$INS_evidence){
          RNA_depth_alt=max(RNA_depth_alt, tmp2$INS)
          RNA_depth_total<-max(RNA_depth_total, tmp2$Depth)
          RNA_evidence<-paste(RNA_evidence, "INS", sep=";")
        }
        RNA_AF=max(0, RNA_depth_alt/RNA_depth_total)
      }
      
    }
    
    
    all_rna<-rbind(all_rna, cbind(RNA_depth_total, RNA_depth_alt, RNA_AF, RNA_evidence))
    
  }
  
  #colnames(all_rna)<-paste0(samp, ".", colnames(all_rna))
  return(all_rna)
}

merge_validation<-function(rna_bam, ref, single.sample.merged, samp){
  rna_val<-get_validation_vars_rna(rna_bam, ref, single.sample.merged, samp)
  
  all_fixed<-single.sample.merged[[4]]
  all_info<-single.sample.merged[[5]]
  all_gt<-single.sample.merged[[6]]
  
  all_variants_discovery<-merge(all_fixed, all_info, by="CHROM_POS_REF_ALT")
  all_variants_discovery<-merge(all_variants_discovery, all_gt, by="CHROM_POS_REF_ALT")
  
  #all_variants_discovery<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]])
  all_filt_variants<-cbind(all_variants_discovery, rna_val)
  
  return(all_filt_variants)
}


