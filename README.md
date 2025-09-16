# capseuratconverter

The `capseuratconverter` package provides functions to convert AnnData (`.h5ad`) files to Seurat v5 objects and saves them as RDS files. The package was created specically for the [Cell Annotation Project](https://celltype.info/) (CAP), so please test it first to see it it meets your requirements. It is written in pure R, with no `reticulate` dependencies. The package currently supports AnnData (`.h5ad`) files created via AnnData v0.10+ and converts the following fields: `X`, `obs`, `var`, `raw.X`, `raw.var`, `obsm`, and `uns`.

## Installation

The package isn't published on CRAN, so it can be installed using `devtools`:

```R
install.packages('devtools')
devtools::install_github('cellannotation/capseuratconverter')
```

Alternatively, you can install it from GitHub releases (the fastest way):

```R
destfile <- "capseuratconverter.tar.gz"
url <- "https://github.com/cellannotation/capseuratconverter/releases/download/v0.5/capseuratconverter_0.5.tar.gz"
download.files(url=url, destfile=destfile)
install.packages(destfile, repos = NULL, type='source')
file.remove(destfile)
```

## Usage

There are two main functions in the package: 
- `h5ad_to_seurat(h5ad_path)` Takes the path to an `.h5ad` file and returns a Seurat v5 object.
- `h5ad2rds(h5ad_path)` Takes the path to an .h5ad file, runs `h5ad_to_seurat`, and saves the object to the same path as the input but replaces the `.h5ad` extension with `.rds`. Returns the path to the RDS file.

## Dependencies

- R version >= 4.0.0
- SeuratObject >= 5.0.0
- Matrix
- rhdf5
- log4r

## Contribution

Fill free to contribute with new issues or PRs! The project is small so there are no special guidlines for that.
