FROM rocker/shiny

EXPOSE 3838

# update package manager & build essentials
RUN apt-get update \
    && apt-get install --yes build-essential

# install dependency required by samtools
RUN apt-get install --yes wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libxml2-dev libcurl4-gnutls-dev libssl-dev libjpeg-dev


# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('shiny', 'shinydashboard','jpeg','Hmisc','BiocManager','githubinstall'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('biomaRt','GenomicFeatures','Gviz','rtracklayer','trackViewer','org.Hs.eg.db','TxDb.Hsapiens.UCSC.hg19.knownGene'))"
RUN R -e "library(githubinstall); githubinstall(c('FusionExpressionPlot'),ask=F)"
RUN R -e "install.packages(c('ggplot2'), repos='http://cran.rstudio.com/')"

ADD app /srv/shiny-server/

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

USER shiny

CMD ["/usr/bin/shiny-server"]