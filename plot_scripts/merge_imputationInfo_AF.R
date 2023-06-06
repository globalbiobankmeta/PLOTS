#!/usr/bin/env Rscript

options(stringsAsFactors=F,bitmapType='cairo')
Sys.setlocale("LC_CTYPE", "C.UTF-8")

library(data.table)
library(optparse)
library(plotrix)
library(ggplot2)

option_list <- list(
        make_option(c("--infile"), type="character", default="",
                help="infile"),
        make_option(c("--AFfile"), type="character", default="",
                help="path to allele frequency file"),
	make_option(c("--imputationInfofile"), type="character", default="",
                help="path to imputation info file"),
	make_option(c("--outfile"), type="character", default="",
                help="path to output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


#opt=list(infile="BioVU_Hirbo_male_VTE_EA_Plink2_10182020.txt.bgz" , AFfile="EA_all_frq.txt", imputationInfofile="biovu_imputationInfo.txt", outfile="BioVU_Hirbo_male_VTE_EA_Plink2_10182020.txt.bgz.EA.txt")

#BioVU_Hirbo_UterineCancer_AA_Plink2_10182020.txt.bgz
#opt=list(infile="BioVU_Hirbo_UterineCancer_AA_Plink2_10182020.txt.bgz" , AFfile="AA_all_frq.txt",imputationInfofile="biovu_imputationInfo.txt", outfile="BioVU_Hirbo_UterineCancer_AA_Plink2_10182020.txt.bgz.AA.txt")

#CHR	BP	ID	REF	ALT	FIRTH?	TEST	OBS_CT	BETA	SE	T_STAT	P
data = fread(paste0("zcat ", opt$infile), header=F, data.table=F)
if(data[1,2] == "BP"){data = data[-1,]}
data$chrpos = paste0(data[,1], ":", data[,2])
colnames(data) = c("#CHR","POS","SNPID","Allele1","Allele2", "FIRTH","TEST","N", "BETA", "SE", "Tstat", "p.value", "chrpos")
#CHR	POS	SNPID	Allele1	Allele2	AC_Allele2	AF_Allele2	imputationInfo	N	BETA	SE	Tstat	p.value	p.value.NA	Is.SPA.converge	varT	varTstar	AF.Cases	AF.Controls	N.Cases	N.Controls
data$SNP2 = paste0(data$SNPID, ":", data$Allele1,":",data$Allele2)
data = data.table(data)
setkey(data,SNP2)

dataAF = fread(opt$AFfile, header=T, data.table=F)
dataAF0 = cbind(dataAF$SNP, dataAF$ALT, dataAF$REF, 1-dataAF$ALT_FREQS)
colnames(dataAF0) = colnames(dataAF)
dataAF2 = rbind(dataAF, dataAF0)
dataAF2$SNP2 = paste0(dataAF2$SNP, ":", dataAF2$REF,":",dataAF2$ALT)
dataAF2 = data.table(dataAF2)
setkey(dataAF2,SNP2)
data1 = merge(data, dataAF2, by.x="SNP2", by.y="SNP2", all.x=T)


dataimpuInfo = fread(opt$imputationInfofile, header=T, data.table=T)
setkey(dataimpuInfo,chrpos)
setkey(data1, chrpos)
data2 = merge(data1, dataimpuInfo, by.x="chrpos", by.y="chrpos", all.x=T)

data2 = data2[which(data2$POS != "BP"),]
data2$AC_Allele2 = as.numeric(data2$ALT_FREQS) * as.numeric(data2$N) *2
data2$AF_Allele2 = data2$ALT_FREQS
data2 = data2[,c("#CHR","POS","SNPID","Allele1","Allele2", "AC_Allele2","AF_Allele2","imputationInfo","N", "BETA", "SE", "Tstat", "p.value")]
write.table(data2, opt$outfile, quote=F, col.names=T, row.names=F, sep="\t")
