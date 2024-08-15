FROM satijalab/seurat:5.0.0

# COPY src/*.R src/
COPY . /app

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

RUN R -e "install.packages('logger')"
RUN R -e "install.packages('usethis')"
RUN R -e "install.packages('devtools')"

WORKDIR /app

CMD ["R"]
