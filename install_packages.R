#!/usr/bin/env Rscript
#install required R packages, from Finnge/SAIGE-IT
req_packages <- c("optparse", "data.table",  "RColorBrewer", "plotrix")
#req_packages <- c("RColorBrewer")
for (pack in req_packages) {
    if(!require(pack, character.only = TRUE)) {
        #install.packages(pack, repos = "https://cloud.r-project.org")
        install.packages(pack, repos='http://cran.us.r-project.org')
    }
}

#devtools::install_github("leeshawn/SKAT")
