#!/usr/bin/env Rscript

options(stringsAsFactors=F,bitmapType='cairo')
Sys.setlocale("LC_CTYPE", "C.UTF-8")

library(data.table)
library(optparse)
library(plotrix)
library(ggplot2)
#opt=list(infile="AF_bothsex_inv_var_meta.gz", bbkpop="EBB_nfe", loofile="./output/leave-EBB_nfe.list")
option_list <- list(
        make_option(c("-b", "--bbkpop"), type="character", default="",
                help="bbk_pop"),
        make_option(c("-m", "--minNingnomAD"), type="numeric", default=100,
                help="100"),
        make_option(c("-i", "--infile"), type="character", default="",
                help="prefix for in sum stats file")
)

        #make_option(c("-o", "--outfile"), type="character", default="",
         #       help="prefix for output file")
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

print(opt)
#opt=list(bbkpop="UKBB_nfe", minNingnomAD=100, infile="UKBB.Karczewski.ThC.ALL.EUR.427.418625.SAIGE.20200802.out.munged.MAC.20.0.INFO.0.3.gz.gnomad_v3_b38_ref_nfe.gz", outfile="UKBB_nfe_gnomAD")

pop=strsplit(opt$bbkpop, split="_")[[1]][2]
print(pop)
##CHR	POS	Allele1	Allele2	AF_Allele2	imputationInfo	BETA	SE	p.value	AF_gnomad_v3_b38_ref_nfe	AF_fc	N
data = fread(paste0("zcat ", opt$infile), header=T, data.table=F, select=c("#CHR", "POS", "Allele1", "Allele2", "AF_Allele2", paste0("AF_gnomad_v3_b38_ref_",pop), "AF_fc", "N"))

#data$notnacount=apply(data[,paste0(loolist[,1], "_AF_Allele2"), drop=F], 1, function(x) sum(!is.na(x)))
#data = data[which(data$notnacount > 1), ]
data = data[which(!is.na(data[,paste0("AF_gnomad_v3_b38_ref_", pop)]) & !is.na(data$AF_Allele2)  ), ]
data$gnomADAF = data[,paste0("AF_gnomad_v3_b38_ref_", pop)]
print(dim(data))
data$AF = data$AF_Allele2
data$AF_new = 1 - data$AF
data$fc = pmax(data$AF,data$gnomADAF)/pmin(data$AF,data$gnomADAF)
data$fc_new = pmax(data$AF_new,data$gnomADAF)/pmin(data$AF_new,data$gnomADAF)

flipIndex = which( ((data$Allele1 == "A" & data$Allele2 == "T") |
        (data$Allele1 == "T" & data$Allele2 == "A") |
        (data$Allele1 == "G" & data$Allele2 == "C") |
        (data$Allele1 == "C" & data$Allele2 == "G") ) & 
(((data$fc > 2) & abs(data$AF - data$gnomADAF) > abs(data$AF_new - data$gnomADAF)) | (data$AF > 0.6 & data$gnomADAF < 0.4) | (data$AF < 0.4 & data$gnomADAF > 0.6) )) 


data$AF[flipIndex] = data$AF_new[flipIndex]
data$fc[flipIndex] = data$fc_new[flipIndex]
data$MAF = pmin(data$AF, 1-data$AF)
print(length(flipIndex))
#rmIndex = which((data$fc > 10 & data$MAF >= 0.01) | (data$fc > 100 & data$MAF < 0.01))
#print(length(rmIndex))
#flipIndex_new = setdiff(flipIndex, rmIndex)
x=data[,c("AF","gnomADAF")]
covM=cov(x)
mdist = mahalanobis(x, center=T, covM)
data$mdist = mdist
data$flip = 0
data$flip[flipIndex] = 1
#flipIndex_norm = setdiff(flipIndex, rmIndex)
write.table(data, paste0(opt$bbkpop,"_postGWASQC_comparetognomAD.txt"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(data[which(data$flipIndex == 1), ], paste0(opt$bbkpop,"_postGWASQC_comparetognomAD_ambiguity.txt"), col.names=T, row.names=F, quote=F, sep="\t")

#data = data[-rmIndex,]
#write.table(data, paste0(opt$bbkpop,"_postGWASQC.txt"), col.names=T, row.names=F, quote=F, sep="\t")
print(dim(data))
title_plot = paste0(nrow(data), " variants in total")
#xl=paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD")
png(paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD_beforeMdistQC.png"), width=1000, height=1000, units="px")
        p <- ggplot(data, aes_string(x="AF", y="gnomADAF")) + xlab(paste0(opt$bbkpop, "_AF_Allele2"))+ylab("AF in gnomAD") + 
          geom_point(alpha=0.1) + 
          theme_minimal(base_size=18) + 
	 ggtitle(title_plot)
        print(p)
dev.off()

png(paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD_squaredMdist.png"), width=1000, height=1000, units="px")
p<-ggplot(data, aes(x=mdist)) + 
  geom_histogram(color="black", fill="white")

p = p + ggtitle(paste0("squared Mahalanobis distance of AF in ", opt$bbkpop, " and gnomAD"))
print(p)
dev.off()

data1 = data[which(data$mdist < 30), ]
print(dim(data1))
data1b = data[which(data$mdist >= 30),]
title_plot = paste0(nrow(data1b), " variants with MD >= 30 were removed")
png(paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD_afterMdistQC_rmMDge30.png"), width=1000, height=1000, units="px")
        p <- ggplot(data1, aes_string(x="AF", y="gnomADAF")) + xlab(paste0(opt$bbkpop, "_AF_Allele2"))+ylab("AF in gnomAD") +
          geom_point(alpha=0.1) +
          theme_minimal(base_size=18) +
	ggtitle(title_plot)
        print(p)
dev.off()
write.table(data1b, paste0(opt$bbkpop,"_postGWASQC_comparetognomAD_MDge30.txt"), col.names=T, row.names=F, quote=F, sep="\t")


data1 = data[order(data$mdist, decreasing=T), ] 
data2 = data1[10001:nrow(data1),]
print(dim(data2))
title_plot = paste0("top 10,000 variants were removed")
png(paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD_afterMdistQC_rmMDtop10k.png"), width=1000, height=1000, units="px")
        p <- ggplot(data2, aes_string(x="AF", y="gnomADAF")) + xlab(paste0(opt$bbkpop, "_AF_Allele2"))+ylab("AF in gnomAD") +
          geom_point(alpha=0.1) +
          theme_minimal(base_size=18) +
        ggtitle(title_plot)
        print(p)
dev.off()
write.table(data1[1:10000,], paste0(opt$bbkpop,"_postGWASQC_comparetognomAD_MDtop10k.txt"), col.names=T, row.names=F, quote=F, sep="\t")

a=floor(nrow(data1)/100)
#data1 = data[order(data$mdist, decreasing=T), ]
data2 = data1[(a+1):nrow(data1),]
print(dim(data2))
title_plot = paste0(a, " variants with top 1 percent MD were removed")
png(paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD_afterMdistQC_rmtop1perc.png"), width=1000, height=1000, units="px")
        p <- ggplot(data2, aes_string(x="AF", y="gnomADAF")) + xlab(paste0(opt$bbkpop, "_AF_Allele2"))+ylab("AF in gnomAD") +
          geom_point(alpha=0.1) +
          theme_minimal(base_size=18) +
        ggtitle(title_plot)
        print(p)
dev.off()
write.table(data1[1:a,], paste0(opt$bbkpop,"_postGWASQC_comparetognomAD_MDtop1perc.txt"), col.names=T, row.names=F, quote=F, sep="\t")

