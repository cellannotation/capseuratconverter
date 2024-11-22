# Development with docker

To build the image run:

```bash
docker build . -t "h5rds"
```

After that it could be run in interactive mode mounting source code `R` directory to support reactive development without rebuilding the docker container:

```bash
docker run -it -v "/your_path_to_R:/app/R" h5rds 
```

Now, one can import the code with:

```R
source("R/h5ad2rds.R")
```

When new changes have been applied, just re-import source code.
