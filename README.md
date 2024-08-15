# capseuratconverter
The dedicated repo for h5ad to Seurat-rds files convertion


# How to run

The docker file is build on `satijalab/seurat:5.0.0`. To buld the image run:

```bash
docker build . -t "h5rds"
```

After that it could be runned in itterative mode mounting directory with code (mounting is needed temporarly to dynamically catch changes in script):

```bash
docker run -it h5rds
```

Or `R` directory could be mounted to support active development without rebuild the docker container:

```bash
docker run -it -v "/your_path_to_R:/R" h5rds 
```

Now one can import the code with R command:

```R
source("src/h5ad2rds.R")
```

