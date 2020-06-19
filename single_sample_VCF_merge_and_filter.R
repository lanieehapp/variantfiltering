#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



##### single sample VCF merger and single sample caller variant filtering 

##INPUT: sample ID and list of file names - HC, S2, DV
##OUTPUT: two sets of variants - one raw merged with filters, one filtered merged with filters

samp<-args[8]
#path<-c("s3://davelab-analysis-results/30427/annovar__haplotypecaller__dna/annovar__haplotypecaller__dna_Annovar.hg38_multianno.vcf", "s3://davelab-analysis-results/30427/annovar__strelka2__dna/annovar__strelka2__dna_Annovar.hg38_multianno.vcf", "s3://davelab-analysis-results/30427/annovar__deepvariant__dna/annovar__deepvariant__dna_Annovar.hg38_multianno.vcf")

path<-c(args[1], args[2], args[3])


source("single_sample_VCF_merge_functions.R")

single.sample.merged<-three.caller.merge(samp, path)



single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])
write.csv(single.sample.all, file=args[4])

single.sample.filt<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]])
write.csv(single.sample.filt, file=args[5])

single.sample.wl<-cbind(single.sample.merged[[7]], single.sample.merged[[8]], single.sample.merged[[9]])
write.csv(single.sample.wl, file=args[6])

save(single.sample.merged, file=args[7])

