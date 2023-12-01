# EMBER creates a unified space for independent breast cancer transcriptomic datasets enabling precision oncology.

This repo contain the scripts used to generate all the figures from the 
EMBER paper and also instructions on how to run the analysis yourself.

## Directory structure

`scripts` contains all the scripts used to generate the images
and associated files. It is organized into a quarto book.

`results` contains all the output of the analysis, including
html and other files. 

`data` contains raw data that is used in the project and is not available
on the server. Note here we are using preprocessed data already available on
the server. For now this folder is empty. 

## Book

The book when generated should be available on the folder `results/book`. 
Also an automatically generated online version is available at [chronchi.github.io/molecular_landscape](https://chronchi.github.io/molecular_landscape).

## Files that are ignored

Overall we are ignoring all the rds, rdb and RData files that are generated
in the analysis. We also ignore the cache folders and the files folders. 
The data folders and results folders are also ignored, they would
make the repo too big. Check the gitignore files for a complete list. The only
folder with its contents that is available from the `results` is the `book`,
which contains the files to visualize the whole analysis.

## Docker

The analysis here has a docker image with all the datasets available. The
docker image is used to automatically run the whole analysis and generate
the report using github actions. Also, if you would like to rerun some
specific parts of the analysis you can download the docker image
`chronchi/molecular_landscape` and run it. This is a Rstudio server
docker image so you will be able to generate the whole analysis in your
computer as well.

### Instructions for running docker

