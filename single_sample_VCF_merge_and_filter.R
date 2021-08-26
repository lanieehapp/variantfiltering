#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



##### single sample VCF merger and single sample caller variant filtering 

##INPUT: sample ID and list of file names - HC, S2, DV
##OUTPUT: three sets of variants - 1) raw merged with filters, 2) filtered merged with filters, 3) whitelist merged 

samp<-args[14]

path<-c(args[1], args[2], args[3])

dna_bam<-args[4]
rna_bam<-args[5]


ref=args[6]
ref_fai=args[7]


source("single_sample_VCF_merge_functions.R")

single.sample.merged<-three.caller.merge(samp, path)

###### return header only if no variants pass filtering #####
if(nrow(single.sample.merged[[4]])==0){
        
        print("No variants pass filtering")
        tmp<-read.csv(text="RNA_depth_total,RNA_depth_alt,RNA_AF,RNA_evidence,Sample_ID")
        all_filt_variants<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]], tmp)
}

if(nrow(single.sample.merged[[4]])>0){
        
        "Variants pass filtering - now running DNA validation"
        
        all_filt_variants<-merge_validation(rna_bam, ref, single.sample.merged, samp)
        
        all_filt_variants <- data.frame(lapply(all_filt_variants, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
        all_filt_variants <- data.frame(lapply(all_filt_variants, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)
        all_filt_variants$Sample_ID<-samp
        
}


all_whitelist<-merge_whitelists(dna_bam, rna_bam, ref, single.sample.merged)


if(nrow(all_whitelist)>0){
        print("Whitelist variants found")
        all_whitelist <- data.frame(lapply(all_whitelist, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
        all_whitelist <- data.frame(lapply(all_whitelist, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)
        all_whitelist$Sample_ID<-samp
        
}





single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])

#put all columns that mess up indels at end of file
bad_cols<-c("HC_BaseQRankSum", "HC_ClippingRankSum", "HC_MQRankSum", "HC_ReadPosRankSum", "S2_SNVHPOL", "S2_CIGAR", "S2_RU", "S2_REFREP", "S2_IDREP", "S2_DPI", "S2_AD", "S2_ADF", "S2_ADR", "S2_FT", "S2_PL", "S2_PS")

#all variants
all_good<-single.sample.all[,!(colnames(single.sample.all) %in% bad_cols)]
all_bad<-single.sample.all[,colnames(single.sample.all) %in% bad_cols]

if("SNVHPOL" %in% colnames(all_bad)){
        colnames(all_bad)[colnames(all_bad)=="SNVHPOL"]<-"S2_SNVHPOL"
        
        correct_order<-cbind(all_bad[,2:5], all_bad[,1], all_bad[,11:15], all_bad[,10], all_bad[,16])
        all_bad<-correct_order
}

single.sample.all<-cbind(all_good, all_bad)

single.sample.all <- data.frame(lapply(single.sample.all, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
single.sample.all <- data.frame(lapply(single.sample.all, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)

single.sample.all$Sample_ID<-samp

#filtered variants
all_good<-all_filt_variants[,!(colnames(all_filt_variants) %in% bad_cols)]
all_bad<-all_filt_variants[,colnames(all_filt_variants) %in% bad_cols]

all_filt_variants<-cbind(all_good, all_bad)

#whitelist
all_good<-all_whitelist[,!(colnames(all_whitelist) %in% bad_cols)]
all_bad<-all_whitelist[,colnames(all_whitelist) %in% bad_cols]

all_whitelist<-cbind(all_good, all_bad)


###filter for exonic, not synonymous, >5 ALT reads
clean_filtered<-all_filt_variants
clean_filtered <- clean_filtered[clean_filtered$Func.refGene=="exonic",]        
clean_filtered <- clean_filtered[!(clean_filtered$ExonicFunc.refGene=="synonymous_SNV"),]

all_depth<-clean_filtered[,grep("*_AD$", colnames(clean_filtered))]
alt_depth<-matrix(NA, nrow=nrow(all_depth), ncol=ncol(all_depth))
for(i in 1:nrow(all_depth)){
        for(j in 1:ncol(all_depth)){
                if(!(is.na(all_depth[i,j]))){
                        alt_depth[i,j]<-unlist(strsplit(all_depth[i,j], ","))[2]
                        
                }
        }
}

alt_depth[is.na(alt_depth)]<-0
max_depth<-apply(alt_depth, 1, max)

clean_filtered<-clean_filtered[max_depth>=5,]


#all variants in .csv
write.table(single.sample.all, file=args[8], row.names  =FALSE, quote=FALSE, sep = "\t")

#filtered variants in .csv
write.table(all_filt_variants, file=args[9], row.names = FALSE, quote = FALSE, sep="\t")

#whitelist variants in .csv
write.table(all_whitelist, file=args[10], row.names=FALSE, quote = FALSE, sep="\t")

write.table(clean_filtered, file=args[15], row.names=FALSE, quote=FALSE, sep="\t")


#simple merging for all variants, filtered variants, whitelist variants
save(single.sample.merged, file=args[11])
#all filtered and whitelist with validation info
#save(all_filt_variants, file="all_filt_variants.RData")
#save(all_whitelist, file="all_whitelist.RData")
save(all_filt_variants, file=args[12])
save(all_whitelist, file=args[13])








