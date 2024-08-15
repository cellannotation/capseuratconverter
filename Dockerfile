FROM r-base:4.4.0

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('rhdf5')"
RUN R -e "BiocManager::install('SeuratObject')"
RUN R -e "BiocManager::install('Matrix')"
RUN R -e "BiocManager::install('logger')"

COPY . /app

WORKDIR /app
