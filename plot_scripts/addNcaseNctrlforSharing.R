#!/usr/bin/env Rscript

options(stringsAsFactors=F)
library(optparse)
library(data.table)

option_list <- list(
  make_option("--metafile", type="character", default="",
    help=""),
  make_option("--bbkfile", type="character", default="",
    help=""),
  make_option("--outfile", type="character", default="",
    help="")
		    
)
## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#infile2=paste0("/net/hunt/disk2/zhowei/project/biobank_meta/meta-analysis-10212020/samplesize/", opt$pheno,"_",opt$sex,".txt.temp")



bbkdata = fread(opt$bbkfile, header=F, data.table=F)
bbkdata = bbkdata[which(bbkdata[,4] > 50),]
bbkpop = paste0(bbkdata[,3], "_", bbkdata[,2])
bbkpopbeta = paste0(bbkpop, "_beta")
#databeta = annodata[,bbkpopbeta]
bbkpopMDgnomAD = paste0(bbkpop, "_MD_gnomad")
#_strandflip
#datapval = annodata[,c(bbkpoppval, "all_inv_var_meta_p")]
selectCol = c("#CHR","POS","REF","ALT","rsid","all_meta_AF", "all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p","all_inv_var_het_p","all_meta_sample_N","all_meta_N", bbkpopMDgnomAD, bbkpopbeta)

annodata = fread(paste0("zcat ", opt$metafile), header=T, data.table=F, sep="\t", select=selectCol)
dataMDgnomAD  = annodata[,bbkpopMDgnomAD]
databeta = annodata[,bbkpopbeta]

countbbk = function(x){
	bbkL = NULL
	nomissingx = names(x)[which(!is.na(x))]
	for(i in nomissingx){
		i2 = strsplit(i, split="_")[[1]][1]
		bbkL = c(bbkL, i2)
	}
	cntbbk = length(unique(bbkL))
	return(cntbbk)
}


cntbbkL = apply(databeta, 1, countbbk)
annodata$Nbbk = cntbbkL


getdirection = function(x){
	dictL = NULL
	for(i in x){
		if(is.na(i)){
			dictL = c(dictL, "?")
		}else{
			if(i >= 0){
				dictL = c(dictL, "+")
			}else{
				dictL = c(dictL, "-")
			}
		}
        }
	dictL2 = paste(dictL, collapse="")
	return(dictL2)
}

direction = apply(databeta, 1,getdirection)
annodata$direction = direction


#strictNovel = function(x){
#	isNovel = FALSE
#	if(x[length(x)] < 5*10^-8){
#		if(sum(x[1:(length(x)-1)] < 5*10^-8, na.rm = T) == 0){
#			isNovel = TRUE
#		}	
#	}
#	return(isNovel)	
#}
#isstrictNovel = apply(datapval, 1, strictNovel)
#annodata$isstrictNovel = isstrictNovel


#countbbk_sig = function(x){
#        bbkL = NULL
#        nomissingx = names(x)[which(!is.na(x))]
#        nomissingxpval = x[which(!is.na(x))]
#        for(i in 1:length(nomissingx)){
#                i2 = strsplit(nomissingx[i], split="_")[[1]][1]
#		if(nomissingxpval[i] < 5*10^-8){
#                	bbkL = c(bbkL, i2)
#		}
#        }
#        cntbbkSig = length(unique(bbkL))
#        return(cntbbkSig)
#}


#cntbbkL_sig = apply(annodata[,c(bbkpoppval)], 1, countbbk_sig)
#annodata$cntbbk_sig = cntbbkL_sig

getQC_MDgnomAD = function(x){
	QC_MDgnomADL = "PASS"
	xnomissing = x[!is.na(x)]
	if(sum(xnomissing == "FAIL") > 0){QC_MDgnomADL="FAIL"}
	return(QC_MDgnomADL)
}	

QC_MDgnomAD = apply(dataMDgnomAD, 1, getQC_MDgnomAD)
annodata$QC_MDgnomAD = QC_MDgnomAD

bbkdata$bbkpop = bbkpop

count_CaseN_CtrlN = function(x, bbkdata_t){
        nomissingx = names(x)[which(!is.na(x))]
        nomissingxName=NULL
        for(i in 1:length(nomissingx)){
        bbk = strsplit(nomissingx[i], split="_")[[1]][1]
        pop = strsplit(nomissingx[i], split="_")[[1]][2]
        nomissingxName = c(nomissingxName, paste0(bbk, "_", pop))
        }


        bbkdata_m = bbkdata_t[which(bbkdata_t$bbkpop %in% nomissingxName), ]
        Ncase = sum(bbkdata_m[,4])
        Nctrl = sum(bbkdata_m[,5])
        return(c(Ncase, Nctrl))
}

CaseN_CtrlN = t(apply(databeta, 1, count_CaseN_CtrlN, bbkdata))

annodata$Ncase = CaseN_CtrlN[,1]
annodata$Nctrl = CaseN_CtrlN[,2]

annodata = annodata[which(annodata$Nbbk > 1), ]
annodata = annodata[,c("#CHR","POS","REF","ALT","rsid","all_meta_AF", "all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p","all_inv_var_het_p","direction",  "Ncase","Nctrl", "all_meta_N", "Nbbk")]
colnames(annodata)[which(colnames(annodata) == "Nbbk")] = "all_meta_N_biobank"
colnames(annodata)[which(colnames(annodata) == "all_meta_N")] = "all_meta_N_dataset"


annodata = annodata[with(annodata, order(annodata[,1], annodata[,2])), ]
write.table(annodata, opt$outfile, sep="\t", quote=F, col.names=T, row.names=F)
