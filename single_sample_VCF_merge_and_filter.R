#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



##### single sample VCF merger and single sample caller variant filtering 

##INPUT: sample ID and list of file names - HC, S2, DV
##OUTPUT: three sets of variants - 1) raw merged with filters, 2) filtered merged with filters, 3) whitelist merged 

samp<-args[9]

path<-c(args[1], args[2], args[3])

dna_bam<-args[4]
rna_bam<-args[5]
ref<-args[6]
ref_fai<-args[7]


whitelist_path<-"filtered_whitelist_08032020.txt"
comm<-paste0('parallel --colsep "\t" samtools mpileup -a -l ', whitelist_path, ' --fasta-ref ',ref, ' ', rna_bam, ' -r {1} :::: ', ref_fai ,' > /data/rna_whitelist.txt'  )

print(comm)
system(comm)

rna_whitelist<-read.csv(file="/data/rna_whitelist.txt", sep="\t", stringsAsFactors = FALSE)


#source("single_sample_VCF_merge_functions.R")

# single.sample.merged<-three.caller.merge(samp, path)
# 
# 
# all_whitelist<-merge_whitelists(dna_bam, rna_bam, ref, single.sample.merged)
# 
# all_filt_variants<-merge_validation(rna_bam, ref, single.sample.merged, samp)
# 
# 
# 
# single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])
# 
# write.csv(single.sample.all, file=args[8])
# 
# #single.sample.filt<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]])
# save(all_filt_variants, file="all_filt_variants.RData")
# write.csv(all_filt_variants, file=args[9])
# 
# save(all_whitelist, file="all_whitelist.RData")
# 
# write.csv(all_whitelist, file=args[10])
# 
# save(single.sample.merged, file=args[10])
# 
# print(args[7])
# 
