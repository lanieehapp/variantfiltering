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

all_whitelist<-merge_whitelists(dna_bam, rna_bam, ref, single.sample.merged)

all_filt_variants<-merge_validation(rna_bam, ref, single.sample.merged, samp)

single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])



#all variants in .csv
write.table(single.sample.all, file=args[8], row.names  =FALSE, quote=FALSE, sep = ",")

#filtered variants in .csv
write.table(all_filt_variants, file=args[9], row.names = FALSE, quote = FALSE, sep=",")

#whitelist variants in .csv
write.table(all_whitelist, file=args[10], row.names=FALSE, quote = FALSE, sep=",")


#simple merging for all variants, filtered variants, whitelist variants
save(single.sample.merged, file=args[11])
#all filtered and whitelist with validation info
#save(all_filt_variants, file="all_filt_variants.RData")
#save(all_whitelist, file="all_whitelist.RData")
save(all_filt_variants, file=args[12])
save(all_whitelist, file=args[13])


