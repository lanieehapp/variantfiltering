##### single sample VCF merger and single sample caller variant filtering 

##INPUT: sample ID and list of file names - HC, S2, DV
##OUTPUT: two sets of variants - one raw merged with filters, one filtered merged with filters

samp<-"21128_T_1_merged"
path<-c("s3://davelab-analysis-results/30427/annovar__haplotypecaller__dna/annovar__haplotypecaller__dna_Annovar.hg38_multianno.vcf", "s3://davelab-analysis-results/30427/annovar__strelka2__dna/annovar__strelka2__dna_Annovar.hg38_multianno.vcf", "s3://davelab-analysis-results/30427/annovar__deepvariant__dna/annovar__deepvariant__dna_Annovar.hg38_multianno.vcf")


source("~/Dropbox (DaveLab)/Lanie_resources/scripts/single_sample_VCF_merge_functions.R")

single.sample.merged<-single.sample.merge(samp, path)



single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])
write.csv(single.sample.all, file="single_sample_all.csv")

single.sample.filt<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]])
write.csv(single.sample.filt, file="single_sample_filt.csv")

single.sample.wl<-cbind(single.sample.merged[[7]], single.sample.merged[[8]], single.sample.merged[[9]])
write.csv(single.sample.wl, file="single_sample_wl.csv")

