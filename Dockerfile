FROM satijalab/seurat:5.0.0

COPY src/*.R src/

RUN R -e "install.packages('logger')"
RUN R -e "install.packages('RestRserve', repos = 'https://cloud.r-project.org')"

CMD ["R"]
