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

# Install capseuratconverter
RUN wget https://github.com/cellannotation/capseuratconverter/releases/download/v0.2/capseuratconverter_0.2.tar.gz -O /tmp/capseuratconverter.tar.gz
RUN R -e "install.packages('/tmp/capseuratconverter.tar.gz', repos = NULL, type='source')"
RUN rm /tmp/capseuratconverter.tar.gz

COPY . /app

WORKDIR /app
