FROM satijalab/seurat:5.0.0

# COPY src/*.R src/
COPY . /app

RUN R -e "install.packages('logger')"
RUN R -e "install.packages('usethis')"

WORKDIR /app

CMD ["R"]
