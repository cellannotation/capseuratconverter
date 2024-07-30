FROM r-base

RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("rhdf5")'

CMD ["R"]
