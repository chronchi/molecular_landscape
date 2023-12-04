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

The analysis here has a docker image with all the datasets available. 
Unfortunately we cannot use github actions to automatically generate
the report and feed to github pages. The docker image size is too big
(~10GB compressed) and the github runners provide up to 14GB SSD storage
space. So instead one would need to run locally the whole docker. For that
before pushing to the main branch of github I check if there is any difference
in the `renv.lock` file. If there is a new image is automatically generated
and submitted to docker hub to get uptodate images for running the analysis. 

After the image is run another Dockerfile is used to generate the report
that will be used in the github pages. The report is saved in the 
docs folder. So if you would like to run the analysis locally the only
thing that you will need to do is run the command below at the root
of this repository. Make sure to change the absolute path to the 
volume so it fetches the scripts as well. 

```bash
# clone the repo to have the latest script available
git clone git@github.com:chronchi/molecular_landscape.git

# build the image that depends on the chronchi/ember:v1 stored
docker build -f Dockerfile.report -t run_analysis .

# Run docker with previously built image to fetch the docs and store in the
# docs folder of the cloned repo
docker run --rm --name generate_docs -v /path/to/github/repo/docs/:/home/rstudio/ember/docs/:rw run_analysis
```

After this you should be able to access the report on `docs/index.html`. 

Moreover, if you want to play with the data and the code, you can access
the Rstudio server available from the docker image directly using the 
commands below. 

``

This will open up the port 8000 (change in the file `` if this port
is being used) and then you can access the rstudio server at
`localhost:8000`.
