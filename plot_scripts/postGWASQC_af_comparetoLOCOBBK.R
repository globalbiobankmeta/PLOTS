#!/usr/bin/env Rscript

options(stringsAsFactors=F,bitmapType='cairo')
Sys.setlocale("LC_CTYPE", "C.UTF-8")

library(data.table)
library(optparse)
library(plotrix)
library(ggplot2)

option_list <- list(
        make_option(c("-b", "--bbkpop"), type="character", default="",
                help="bbk_pop"),
        make_option(c("-l", "--loofile"), type="character", default="",
                help="prefix for output files"),
        make_option(c("-i", "--infile"), type="character", default="",
                help="prefix for in sum stats file"),
        make_option(c("-o", "--outfile"), type="character", default="",
                help="prefix for output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

print(opt)

loolist = fread(opt$loofile, header=F, data.table=F)
print(loolist[,1])

data = fread(paste0("zcat ", opt$infile), header=T, data.table=F, select=c("#CHR", "POS", "REF", "ALT", paste0(opt$bbkpop, "_AF_Allele2"), paste0(loolist[,1], "_AF_Allele2")))

data$notnacount=apply(data[,paste0(loolist[,1], "_AF_Allele2"), drop=F], 1, function(x) sum(!is.na(x)))
data = data[which(data$notnacount > 1), ]

data = data[which(!is.na(data[,paste0(opt$bbkpop, "_AF_Allele2")])), ]
print(dim(data))
data$otherBiobankFreq = rowMeans(data[,paste0(loolist[,1], "_AF_Allele2"),drop=F], na.rm=T)
data$AF = data[,paste0(opt$bbkpop, "_AF_Allele2")]
data$fc = pmax(data$AF,data$otherBiobankFreq)/(pmin(data$AF,data$otherBiobankFreq))
data$AF_new = 1 - data[,paste0(opt$bbkpop, "_AF_Allele2")]
data$fc_new = pmax(data$AF_new,data$otherBiobankFreq)/(pmin(data$AF_new,data$otherBiobankFreq))

flipIndex = which( ((data$REF == "A" & data$ALT == "T") |
        (data$REF == "T" & data$ALT == "A") |
        (data$REF == "G" & data$ALT == "C") |
        (data$REF == "C" & data$ALT == "G") ) & 
(((data$fc > 2) & abs(data$AF - data$otherBiobankFreq) > abs(data$AF_new - data$otherBiobankFreq)) | (data$AF > 0.6 & data$otherBiobankFreq < 0.4) | (data$AF < 0.4 & data$otherBiobankFreq > 0.6) )) 


data$AF[flipIndex] = data$AF_new[flipIndex]
data$fc[flipIndex] = data$fc_new[flipIndex]
data$MAF = pmin(data$AF, 1-data$AF)
print(length(flipIndex))
rmIndex = which((data$fc > 10 & data$MAF >= 0.01) | (data$fc > 100 & data$MAF < 0.01))
print(length(rmIndex))

flipIndex_new = setdiff(flipIndex, rmIndex)

#flipIndex_norm = setdiff(flipIndex, rmIndex)
write.table(data[flipIndex_new, c("#CHR", "POS", "REF", "ALT", paste0(opt$bbkpop, "_AF_Allele2"), "otherBiobankFreq")], paste0(opt$bbkpop,"_postGWASQC_ambiguity.txt"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(data[rmIndex, c("#CHR", "POS", "REF", "ALT", paste0(opt$bbkpop, "_AF_Allele2"), "otherBiobankFreq", "fc")], paste0(opt$bbkpop,"_postGWASQC_remove.txt"), col.names=T, row.names=F, quote=F, sep="\t")
data = data[-rmIndex,]
#write.table(data, paste0(opt$bbkpop,"_postGWASQC.txt"), col.names=T, row.names=F, quote=F, sep="\t")
print(dim(data))
xl=paste0(opt$bbkpop,"_postGWASQCd")
png(paste0(opt$bbkpop,"_postGWASQCd.png"), width=1000, height=1000, units="px")
        p <- ggplot(data, aes_string(x="AF", y="otherBiobankFreq")) + xlab(paste0(opt$bbkpop, "_AF_Allele2"))+ylab("average AF in other biobanks") + 
          geom_point(alpha=0.1) + 
          theme_minimal(base_size=18)
        print(p)
dev.off()

