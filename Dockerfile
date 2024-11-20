FROM r-base:4.4.0

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('rhdf5')"
RUN R -e "BiocManager::install('SeuratObject')"
RUN R -e "BiocManager::install('Matrix')"
RUN R -e "BiocManager::install('log4r')"


RUN apt-get update && apt-get install -y \
    libfontconfig1-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    build-essential \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('cellannotation/capseuratconverter')"

COPY . /app

WORKDIR /app
