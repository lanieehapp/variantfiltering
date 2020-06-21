##single sample VCF merge functions


library(vcfR)
library(stringr)
library(plyr)
library(parallel)
library(doParallel)
library(data.table)
library(qqman)
library(bedr)


single.sample.merge<-function(samp, path){

  
  file.list<-c("~/tmp/tmp_HC.vcf.gz", "~/tmp/tmp_S2.vcf.gz", "~/tmp/tmp_DV.vcf.gz")
  
  single.sample.merged<-three.caller.merge(samp, file.list)
  
  return(single.sample.merged)
}

get_vcf = function(file=file, caller=caller) {
  if(startsWith(file, "gs")){
    comm = paste("/Users/laniehapp/google-cloud-sdk/bin/gsutil cp ",file, " ~/tmp/tmp_", caller,".vcf.gz",sep="")
    
  }
  if(startsWith(file, "s3")){
    comm = paste("/usr/local/bin/aws s3 cp ", file, " ~/tmp/tmp_", caller,".vcf.gz",sep="")
  }
  print(comm)
  system(comm)
}

three.caller.merge<-function(samp, file.list){
  
  all.data<-read.vcfR(file=file.list[1])
  all.hc.data<-parse.hc.vcf(all.data)
  
  all.data<-read.vcfR(file=file.list[2])
  all.s2.data<-parse.s2.vcf(all.data)
  
  all.data<-read.vcfR(file=file.list[3])
  all.dv.data<-parse.dv.vcf(all.data)
  
  all.data.merged<-merge.vcfs(samp,all.hc.data, all.dv.data, all.s2.data)
  gt.summary<-get.gt.summary(all.data.merged)
  all.fix.merged<-all.data.merged[[1]]
  all.info.merged<-all.data.merged[[2]]
  all.gt.merged<-all.data.merged[[3]]
  
  all.gt.merged<-data.frame(all.gt.merged[,1],gt.summary, all.gt.merged[,2:ncol(all.gt.merged)], stringsAsFactors = FALSE, check.names = FALSE)
  colnames(all.gt.merged)[1]<-"CHROM_POS_REF_ALT"
  
  
  #put fix in same order as gt and info
  all.fix.merged<-all.fix.merged[order(all.fix.merged$CHROM_POS_REF_ALT),]
  
  ##remove deep variant double entry variants
  if(nrow(all.fix.merged) != nrow(all.gt.merged)){
    tmp<-data.frame(table(all.gt.merged$CHROM_POS_REF_ALT))
    bad.vars<-as.character(tmp[tmp$Freq>1,1])
  }
  if(exists("bad.vars")){
    for( var in bad.vars){
      all.info.merged<-all.info.merged[!grepl(var, all.info.merged$CHROM_POS_REF_ALT),]
      all.gt.merged<-all.gt.merged[!grepl(var, all.gt.merged$CHROM_POS_REF_ALT),]
      all.fix.merged<-all.fix.merged[!grepl(var, all.fix.merged$CHROM_POS_REF_ALT),]
    }
  }
  
  all.info.merged<-all.info.merged[!grepl("ERCC", all.info.merged$CHROM_POS_REF_ALT),]
  all.gt.merged<-all.gt.merged[!grepl("ERCC", all.gt.merged$CHROM_POS_REF_ALT),]
  all.fix.merged<-all.fix.merged[!grepl("ERCC", all.fix.merged$CHROM),]
  
  all.info.merged<-all.info.merged[!grepl("GL000251.2", all.info.merged$CHROM_POS_REF_ALT),]
  all.gt.merged<-all.gt.merged[!grepl("GL000251.2", all.gt.merged$CHROM_POS_REF_ALT),]
  all.fix.merged<-all.fix.merged[!grepl("GL000251.2", all.fix.merged$CHROM),]
  
  all.info.merged<-cbind(all.info.merged, get.annovar.filters(all.fix.merged))
  

  
  ##apply basic filter and caller filter
  all.fix.filt<-all.fix.merged[all.info.merged$Caller_Filter & all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
  all.info.filt<-all.info.merged[all.info.merged$Caller_Filter & all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
  all.gt.filt<-all.gt.merged[all.info.merged$Caller_Filter & all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
  
  
  ##calculate filter for good, rare variants in windows 
  all.vars<-NULL
  for(chr in unique(all.fix.filt$CHROM)){
    
    vars<-all.fix.filt[all.fix.filt$CHROM==chr,]
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
  
  all.vars<-all.vars[order(all.vars$CHROM_POS_REF_ALT),]
  all.info.filt$CLUST<-all.vars$CLUST
  #colnames(all.info.merged)[ncol(all.info.merged)]<-paste0(samp, ".CLUST")
  
  ##AF filter
  all.info.filt$AF_Filter<-all.gt.filt$afMax<0.25
  
  ##depth filter
  all.info.filt$Depth_Filter<-all.gt.filt$dpMax>5
  
  
  all.fix.filt<-all.fix.filt[!(all.info.filt$CLUST & all.info.filt$AF_Filter) & all.info.filt$Depth_Filter,]
  all.gt.filt<-all.gt.filt[!(all.info.filt$CLUST & all.info.filt$AF_Filter) & all.info.filt$Depth_Filter,]
  all.info.filt<-all.info.filt[!(all.info.filt$CLUST & all.info.filt$AF_Filter) & all.info.filt$Depth_Filter,]
  
  all.fix.wl<-all.fix.merged[all.info.merged$whitelist,]
  all.info.wl<-all.info.merged[all.info.merged$whitelist,]
  all.gt.wl<-all.gt.merged[all.info.merged$whitelist,]
  
  
  
  colnames(all.gt.merged)[2:ncol(all.gt.merged)]<-paste(samp, colnames(all.gt.merged)[2:ncol(all.gt.merged)], sep=".")
  colnames(all.info.merged)[2:ncol(all.info.merged)]<-paste(samp, colnames(all.info.merged)[2:ncol(all.info.merged)], sep=".")
  
  colnames(all.gt.filt)[2:ncol(all.gt.filt)]<-paste(samp, colnames(all.gt.filt)[2:ncol(all.gt.filt)], sep=".")
  colnames(all.info.filt)[2:ncol(all.info.filt)]<-paste(samp, colnames(all.info.filt)[2:ncol(all.info.filt)], sep=".")
  
  colnames(all.gt.wl)[2:ncol(all.gt.wl)]<-paste(samp, colnames(all.gt.wl)[2:ncol(all.gt.wl)], sep=".")
  colnames(all.info.wl)[2:ncol(all.info.wl)]<-paste(samp, colnames(all.info.wl)[2:ncol(all.info.wl)], sep=".")
  
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
  
  all.fix.hc$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  all.info.hc$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  all.gt.parsed.hc$CHROM_POS_REF_ALT<-CHROM_POS_REF_ALT
  
  
  
  
  
  
  
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
  
  
  
  samp.meta<-data.frame(ncallers,unlist(af.max), dp.max)
  
  colnames(samp.meta)<-c("nCallers", "afMax", "dpMax")
  return(samp.meta)
  
}

get.annovar.filters<-function(all.fix.merged){
  pop.freq.data<-suppressWarnings(sapply(all.fix.merged[,grep("gnomAD|ExAC|esp6500",colnames(all.fix.merged),value=T)], as.numeric))
  pop.freq.data[is.na(pop.freq.data)] <- 0
  all.fix.merged$pop.freq.max.all <- apply(pop.freq.data,1,max)
  rm(pop.freq.data)
  
  all.filters<-data.frame(all.fix.merged$CHROM_POS_REF_ALT, stringsAsFactors = FALSE)
  all.filters$exonic <- all.fix.merged$Func.refGene=="exonic"
  load("whitelist.RData")
  all.filters$whitelist<-all.fix.merged$CHROM_POS_REF_ALT %in% whitelist$CHROM_POS_REF_ALT
  all.filters$rare.variants <-all.fix.merged$pop.freq.max.all<=0.001
  all.filters$not.syn <- all.fix.merged[,'ExonicFunc.refGene']!="synonymous_SNV" 
  all.filters$not.superdups <- all.fix.merged$genomicSuperDups=="."
  all.filters$basic.filters<- (all.filters$whitelist) | (all.filters$rare & all.filters$not.syn & all.filters$not.superdups)  

  
  
  ###repeatMasker
 

  load("repeat_masker.RData")
  options(scipen = 999)
  var.locs<-c(paste0(all.fix.merged[,1], ":", all.fix.merged[,2],"-", as.numeric(all.fix.merged[,2])+1))
  #var.locs<-cbind(all.fix.merged[,1:2], all.fix.merged[,2]+1)
  table(is.valid.region(var.locs))
  var.locs<-bedr.sort.region(var.locs)
  
  not.repeatmasker<-!(in.region(var.locs, bed.file.filt))
  
  names(not.repeatmasker)<-var.locs
  not.repeatmasker<-not.repeatmasker[order(names(not.repeatmasker))]
  all.filters$not.repeatmasker<-not.repeatmasker

  
  return(all.filters[,2:8])
  
}
