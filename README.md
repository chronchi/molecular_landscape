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

## Docs

The analysis when generated is available on the folder `docs`. 
Also an online version is available at [chronchi.github.io/molecular_landscape](https://chronchi.github.io/molecular_landscape).

## Files that are ignored

Overall we are ignoring all the rds, rdb and RData files that are generated
in the analysis. We also ignore the cache folders and the files folders. 
The data folders and results folders are also ignored, they would
make the repo too big. Check the gitignore files for a complete list. 

## Docker

The analysis here has a docker image with all the datasets available. 
Unfortunately we cannot use github actions to automatically generate
the report and feed to github pages. The docker image size is too big
(~10GB compressed) and the github runners provide up to 14GB SSD storage
space. So instead one would need to run locally the whole docker. For that
before pushing to the main branch of github I check if there is any difference
in the `renv.lock` file. If there is, a new image is automatically generated
and submitted to docker hub to provide images for running the analysis. 

One step that was crucial when creating the image was to first 
isolate the cache from the package directly into the
private library (`renv::isolate()`). This way we ensure that all
the packages are copied to the docker image.

After the image is run another Dockerfile is used to generate the report
that will be used in the github pages. The report is saved in the 
docs folder. So if you would like to run the analysis locally the only
thing that you will need to do is run the command below at the root
of this repository. You will need sudo access to run it.

```bash
# clone the repo to have the latest script available
git clone git@github.com:chronchi/molecular_landscape.git

bash generate_docs.sh
```

After this you should be able to access the report on `docs/index.html`. 

The docker image does not contain all the intermediate files necessary
to run the analysis. They are generated when creating the docs.

Moreover, if you want to play with the data and the code, you can access
the RStudio server available from the docker image directly using the 
commands below. The username is `rstudio` and the password is 
`ember`. The scripts in the docker image are the latest available upon
the creation of the image. Moreover, the rds files associated with
each chapter are also available. 

```bash
docker run \
    -p 8000:8787 \
    --name ember \
    -e PASSWORD=ember \
    chronchi/ember:latest
```

RStudio can be run from `localhost:8000` after that.
