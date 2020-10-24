FROM rocker/r-ver:3.6.1

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget


# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    chown shiny:shiny /var/lib/shiny-server

EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]

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
