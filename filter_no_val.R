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

#all_whitelist<-merge_whitelists(dna_bam, rna_bam, ref, single.sample.merged)

#all_filt_variants<-merge_validation(rna_bam, ref, single.sample.merged, samp)

single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])
all_whitelist<-cbind(single.sample.merged[[7]], single.sample.merged[[8]], single.sample.merged[[9]])
all_filt_variants<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]])





single.sample.all <- data.frame(lapply(single.sample.all, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
single.sample.all <- data.frame(lapply(single.sample.all, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)

all_whitelist <- data.frame(lapply(all_whitelist, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
all_whitelist <- data.frame(lapply(all_whitelist, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)

all_filt_variants <- data.frame(lapply(all_filt_variants, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
all_filt_variants <- data.frame(lapply(all_filt_variants, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)


######put all variants, filtered variants, and whitelist variants into MAF format #########3

###all variants###

# n<-nrow(single.sample.all)
# single.sample.all$REF<-as.character(single.sample.all$REF)
# single.sample.all$ALT<-as.character(single.sample.all$ALT)
# end_pos<-unlist(lapply(1:n, function(x){as.numeric(single.sample.all$POS[x]) + max(nchar(single.sample.all$REF[x]), nchar(single.sample.all$ALT[x]))}))-1 
# 
# 
# single.sample.all$AAChange.refGene<-as.character(single.sample.all$AAChange.refGene)
# aa_change<-strsplit(single.sample.all$AAChange.refGene, ":")
# c_change<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[2], x[4], sep=":"), paste("none"))}))
# p_change<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[2], x[5], sep=":"), paste("none"))}))
# exon_num<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[3]), paste("none"))}))
# 
# gnomad_genome<-single.sample.all[,grep("gnomad_genome_", colnames(single.sample.all))]
# gnomad_genome<-apply(gnomad_genome, 2, function(x){as.numeric(paste(x))})
# gnomad_genome_max<-apply(gnomad_genome, 1, max, na.rm=TRUE)
# gnomad_genome_max[gnomad_genome_max==-Inf]<-0
# var.type<-rep(NA, n)
# var.type[nchar(single.sample.all$REF)>1]<-"DEL"
# var.type[nchar(single.sample.all$ALT)>1]<-"INS"
# var.type[is.na(var.type)]<-"SNP"
# 
# all.maf<-data.frame(single.sample.all$Gene.refGene, rep("none", n), rep("davelab", n), rep("hg38", n), single.sample.all$CHROM,
#                     single.sample.all$POS, end_pos, rep("+", n), single.sample.all$Func.refGene, single.sample.all$ExonicFunc.refGene, var.type, 
#                     single.sample.all$REF,  single.sample.all$ALT, rep("none", n), single.sample.all$avsnp150, rep(samp, n), c_change, 
#                     p_change, exon_num, single.sample.all[,grep(".nCallers", colnames(single.sample.all))],single.sample.all[,grep(".dpMax", colnames(single.sample.all))], 
#                     single.sample.all[,grep(".afMax", colnames(single.sample.all))], single.sample.all$gnomad_exome_AF, single.sample.all$gnomad_exome_AF_popmax, 
#                     single.sample.all$gnomad_genome_AF, gnomad_genome_max, single.sample.all[, grep("pop.freq.max.all", colnames(single.sample.all))], single.sample.all$CADD_phred,
#                     single.sample.all$cadd16gt10, single.sample.all$Interpro_domain, single.sample.all$CLNSIG, stringsAsFactors = FALSE)
# 
# colnames(all.maf)<-c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", 
#                      "Variant_Location", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
#                      "dbSNP_RS", "Tumor_Sample_ID", "HVGSc", "HVGSp", "Exon_Number", "n_callers", "t_max_caller_depth", "t_max_caller_af", "gnomad_exome_all", 
#                      "gnomad_exome_popmax", "gnomad_genome_all", "gnomad_genome_popmax", "all_db_popmax", "cadd_phred", "cadd16_gt10", "interpro_domain", 
#                      "clinical_significance")
# 
# ###filtered variants###
# n<-nrow(all_filt_variants)
# 
# all_filt_variants$REF<-as.character(all_filt_variants$REF)
# all_filt_variants$ALT<-as.character(all_filt_variants$ALT)
# 
# end_pos<-unlist(lapply(1:n, function(x){as.numeric(all_filt_variants$POS[x]) + max(nchar(all_filt_variants$REF[x]), nchar(all_filt_variants$ALT[x]))}))-1 
# 
# all_filt_variants$AAChange.refGene<-as.character(all_filt_variants$AAChange.refGene)
# aa_change<-strsplit(all_filt_variants$AAChange.refGene, ":")
# c_change<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[2], x[4], sep=":"), paste("none"))}))
# p_change<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[2], x[5], sep=":"), paste("none"))}))
# exon_num<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[3]), paste("none"))}))
# 
# gnomad_genome<-all_filt_variants[,grep("gnomad_genome_", colnames(all_filt_variants))]
# gnomad_genome<-apply(gnomad_genome, 2, function(x){as.numeric(paste(x))})
# gnomad_genome_max<-apply(gnomad_genome, 1, max, na.rm=TRUE)
# gnomad_genome_max[gnomad_genome_max==-Inf]<-0
# var.type<-rep(NA, n)
# var.type[nchar(all_filt_variants$REF)>1]<-"DEL"
# var.type[nchar(all_filt_variants$ALT)>1]<-"INS"
# var.type[is.na(var.type)]<-"SNP"
# 
# filt.maf<-data.frame(all_filt_variants$Gene.refGene, rep("none", n), rep("davelab", n), rep("hg38", n), all_filt_variants$CHROM,
#                      all_filt_variants$POS, end_pos, rep("+", n), all_filt_variants$Func.refGene, all_filt_variants$ExonicFunc.refGene, var.type, 
#                      all_filt_variants$REF,  all_filt_variants$ALT, rep("none", n), all_filt_variants$avsnp150, rep(samp, n), c_change, 
#                      p_change, exon_num, all_filt_variants[,grep(".nCallers", colnames(all_filt_variants))],all_filt_variants[,grep(".dpMax", colnames(all_filt_variants))], 
#                      all_filt_variants[,grep(".afMax", colnames(all_filt_variants))], all_filt_variants$gnomad_exome_AF, all_filt_variants$gnomad_exome_AF_popmax, 
#                      all_filt_variants$gnomad_genome_AF, gnomad_genome_max, all_filt_variants[, grep("pop.freq.max.all", colnames(all_filt_variants))], all_filt_variants$CADD_phred,
#                      all_filt_variants$cadd16gt10, all_filt_variants$Interpro_domain, all_filt_variants$CLNSIG, all_filt_variants[,grep(".RNA_", colnames(all_filt_variants))], stringsAsFactors = FALSE)
# 
# colnames(filt.maf)<-c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", 
#                       "Variant_Location", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
#                       "dbSNP_RS", "Tumor_Sample_ID", "HVGSc", "HVGSp", "Exon_Number", "n_callers", "t_max_caller_depth", "t_max_caller_af", "gnomad_exome_all", 
#                       "gnomad_exome_popmax", "gnomad_genome_all", "gnomad_genome_popmax", "all_db_popmax", "cadd_phred", "cadd16_gt10", "interpro_domain", 
#                       "clinical_significance", "validation_RNA_total_depth", "validation_RNA_alt_depth", "validation_RNA_af", "validation_RNA_evidence")
# 
# 
# # ###whitelist variants###
# n<-nrow(all_whitelist)
# 
# #end_pos<-unlist(lapply(1:n, function(x){as.numeric(all_whitelist$Start[x]) + max(nchar(all_whitelist$Ref[x]), nchar(all_whitelist$Alt[x]))}))-1
# 
# all_whitelist$AAChange.refGene.x<-as.character(all_whitelist$AAChange.refGene.x)
# aa_change<-strsplit(all_whitelist$AAChange.refGene.x, ":")
# c_change<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[2], x[4], sep=":"), paste("none"))}))
# p_change<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[2], x[5], sep=":"), paste("none"))}))
# exon_num<-unlist(lapply(aa_change, function(x){ifelse(length(x)>1, paste(x[3]), paste("none"))}))
# 
# gnomad_genome<-all_whitelist[,grep("gnomad_genome_", colnames(all_whitelist))]
# gnomad_genome<-apply(gnomad_genome, 2, function(x){as.numeric(paste(x))})
# gnomad_genome_max<-apply(gnomad_genome, 1, max, na.rm=TRUE)
# gnomad_genome_max[gnomad_genome_max==-Inf]<-0
# 
# all_whitelist$Ref<-as.character(all_whitelist$Ref)
# all_whitelist$Alt<-as.character(all_whitelist$Alt)
# var.type<-rep(NA, n)
# var.type[nchar(all_whitelist$Ref)>1]<-"DEL"
# var.type[nchar(all_whitelist$Alt)>1]<-"INS"
# var.type[is.na(var.type)]<-"SNP"
# 
# wl.maf<-data.frame(all_whitelist$Gene.refGene.x, rep("none", n), rep("davelab", n), rep("hg38", n), all_whitelist$Chr,
#                    all_whitelist$Start, all_whitelist$End, rep("+", n), all_whitelist$Func.refGene.x, all_whitelist$ExonicFunc.refGene.x, var.type,
#                    all_whitelist$Ref,  all_whitelist$Alt, rep("none", n), all_whitelist$avsnp150.x, rep(samp, n), c_change,
#                    p_change, exon_num, all_whitelist[,grep(".nCallers", colnames(all_whitelist))],all_whitelist[,grep(".dpMax", colnames(all_whitelist))],
#                    all_whitelist[,grep(".afMax", colnames(all_whitelist))], all_whitelist$gnomad_exome_AF.x, all_whitelist$gnomad_exome_AF_popmax,
#                    all_whitelist$gnomad_genome_AF.x, gnomad_genome_max, all_whitelist[, grep("pop.freq.max.all", colnames(all_whitelist))], all_whitelist$CADD_phred,
#                    all_whitelist$cadd16gt10, all_whitelist$Interpro_domain.x, all_whitelist$CLNSIG.x, all_whitelist[,grep(".DNA_", colnames(all_whitelist))], all_whitelist[,grep(".RNA_", colnames(all_whitelist))], stringsAsFactors = FALSE)
# 
# colnames(wl.maf)<-c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand",
#                     "Variant_Location", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
#                     "dbSNP_RS", "Tumor_Sample_ID", "HVGSc", "HVGSp", "Exon_Number", "n_callers", "t_max_caller_depth", "t_max_caller_af", "gnomad_exome_all",
#                     "gnomad_exome_popmax", "gnomad_genome_all", "gnomad_genome_popmax", "all_db_popmax", "cadd_phred", "cadd16_gt10", "interpro_domain",
#                     "clinical_significance","validation_DNA_total_depth", "validation_DNA_alt_depth", "validation_DNA_af", "validation_RNA_total_depth", "validation_RNA_alt_depth", "validation_RNA_af")
# 



####################################################################
###### results output ##############################################
####################################################################


#all variants in .csv
write.table(single.sample.all, file=args[8], row.names  =FALSE, quote=FALSE, sep = "\t")

#filtered variants in .csv
write.table(all_filt_variants, file=args[9], row.names = FALSE, quote = FALSE, sep="\t")

#whitelist variants in .csv
write.table(all_whitelist, file=args[10], row.names=FALSE, quote = FALSE, sep="\t")


#simple merging for all variants, filtered variants, whitelist variants
save(single.sample.merged, file=args[11])
#all filtered and whitelist with validation info
#save(all_filt_variants, file="all_filt_variants.RData")
#save(all_whitelist, file="all_whitelist.RData")
save(all_filt_variants, file=args[12])
save(all_whitelist, file=args[13])

#maf files
# write.table(all.maf, file=args[15], row.names=FALSE, quote=FALSE, sep="\t")
# write.table(filt.maf, file=args[16], row.names=FALSE, quote=FALSE, sep="\t")
# write.table(wl.maf, file=args[17], row.names=FALSE, quote=FALSE, sep="\t")
# 
# system("cp dna_whitelist.txt /data/")
# system("cp rna_whitelist.txt /data/")
# system("cp rna_filt.txt /data/")







