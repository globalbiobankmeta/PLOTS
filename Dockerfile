FROM gcr.io/finngen-refinery-dev/bioinformatics:0.6

RUN apt-get install -y libjpeg-dev

ADD install_packages.R /usr/local/bin/
RUN Rscript /usr/local/bin/install_packages.R

ADD plot_scripts /plot_scripts
RUN Rscript /plot_scripts/survPlot.r --help
RUN Rscript /plot_scripts/QQplot.r --help
RUN Rscript /plot_scripts/ManhattanPlot.r --help
RUN /plot_scripts/postGWASQC_af_comparetoLOCOBBK.R --help
RUN /plot_scripts/mergeImputationInfo.R --help

RUN /plot_scripts/postGWASQC_af_compareto_gnomAD.r --help


RUN /plot_scripts/mergePhenoGeno.r --help
RUN /plot_scripts/survPlot.r --help

RUN /plot_scripts/evaluateMD_compareto_gnomAD.r --help
RUN /plot_scripts/addNcaseNctrlforSharing.R --help
RUN /plot_scripts/survPlot.r --help
