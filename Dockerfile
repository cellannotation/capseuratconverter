FROM satijalab/seurat:5.0.0

COPY src/*.R src/

RUN R -e "install.packages('logger')"

CMD ["R"]

# FROM r-base

# RUN R -e 'install.packages("BiocManager")'
# RUN R -e 'BiocManager::install("rhdf5")'

