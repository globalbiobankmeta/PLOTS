#!/usr/bin/env Rscript


options(stringsAsFactors=F,bitmapType='cairo')
Sys.setlocale("LC_CTYPE", "C.UTF-8")

library(data.table)
library(optparse)
library(plotrix)
library(ggplot2)

option_list <- list(
        make_option(c("-i", "--infofile"), type="character", default="",
                help="prefix for info files"),
        make_option(c("-f", "--infile"), type="character", default="",
                help="prefix for in sum stats file"),
        make_option(c("-o", "--outfile"), type="character", default="",
                help="prefix for output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

print(opt)

infodata = fread(opt$infofile, header=T, data.table=T, key="ID", select=c("ID","R2"))
print(head(infodata))


data = fread(paste0("zcat ", opt$infile), header=T, data.table=T, key="SNPID")
print(head(data))

datanew = merge(data, infodata, by.x="SNPID", by.y="ID",all.x=T)

setnames(datanew,"R2","imputationInfo")


fwrite(datanew, opt$outfile, quote=F, row.names=F, col.names=T, sep="\t")

