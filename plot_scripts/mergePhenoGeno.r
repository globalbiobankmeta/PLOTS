#!/usr/bin/env Rscript
#options("gwasurvivr.cores"=8)
library(data.table)
require(optparse) #install.packages("optparse")
print(sessionInfo())

option_list <- list(
  make_option("--pheno", type="character", default="",
    help="pheno"),
  make_option("--phenofile", type="character", default="",
    help="pheno file"),
  make_option("--infile", type="character", default="",
    help="infile")
)
  #make_option("--outfile", type="character", default="",
  #  help="path to output file")

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#phenofolder="/home/wei/survival/realData/FinnGen/pheno/output_r5/withBirthYear/withBirthYear/"
#phenofile=paste0("r6_", opt$pheno,".pheno.txt")
phenodata = fread(opt$phenofile, header=T)
phenodata = data.frame(phenodata)
#E4_DM2_STRICT_H7_RETINOPATHYDIAB_survival_p_lt_5e8.txt.markerlist.all.raw
#infile=paste0("/home/wei/survival/realData/FinnGen/output_PoissonGLMM/download/merged/",opt$pheno,"_survival_p_lt_5e8.txt.markerlist.all.raw")
indata = fread(opt$infile, header=T)
indata = data.frame(indata)
for (i in 7:ncol(indata)){
	a=colnames(indata)[i]
	b=strsplit(a, split="_")[[1]]
	b2 = paste0(b[1],"_",b[2],"_",b[3],"_",b[4])
	colnames(indata)[i] = b2
}

datamerge = merge(phenodata, indata, by.x="FINNGENID", by.y="IID")
write.table(datamerge, paste0(opt$infile, ".",opt$pheno,".txt"), quote=F, col.names=T, row.names=F)


