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
data = fread(paste0("zcat ", opt$infile), header=T, data.table=F, select=c("#CHR", "POS", "Allele1", "Allele2", paste0("AF_gnomad_v3_b38_ref_",pop), "AF", "N", "mdist"))

#data$notnacount=apply(data[,paste0(loolist[,1], "_AF_Allele2"), drop=F], 1, function(x) sum(!is.na(x)))
#data = data[which(data$notnacount > 1), ]
#data = data[which(!is.na(data[,paste0("AF_gnomad_v3_b38_ref_", pop)]) & !is.na(data$AF_Allele2)  ), ]
data$gnomADAF = data[,paste0("AF_gnomad_v3_b38_ref_", pop)]
print(dim(data))

#a = c(opt$bbkpop, mean(data$mdist), sd(data$mdist), mean(data$mdist)+3*(sd(data$mdist)))
#write.table(t(a), paste0(opt$bbkpop,"_squaremd_summary.txt"), col.names=T, row.names=F, quote=F, sep="\t")


cutoff=mean(data$mdist)+3*(sd(data$mdist))
data1 = data[which(data$mdist < cutoff), ]
print(dim(data1))
data1b = data[which(data$mdist >= cutoff),]
#title_plot = paste0(nrow(data1b), " variants with MD >= mean+3*SD ", cutoff, " were removed")
#png(paste0(opt$bbkpop,"_postGWASQCd_comparetognomAD_afterMdistQC_rmMDgt3SDfromMean.png"), width=1000, height=1000, units="px")
#        p <- ggplot(data1, aes_string(x="AF", y="gnomADAF")) + xlab(paste0(opt$bbkpop, "_AF_Allele2"))+ylab("AF in gnomAD") +
#          geom_point(alpha=0.1) +
#          theme_minimal(base_size=18) +
#	ggtitle(title_plot)
#        print(p)
#dev.off()
write.table(data1b, paste0(opt$bbkpop,"_postGWASQC_comparetognomAD_MDge3SDfromMean.txt"), col.names=T, row.names=F, quote=F, sep="\t")
