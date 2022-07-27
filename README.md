# Molecular landscape: a new framework to think about personalized medicine

This repo contain the scripts used to generate the images from the 
molecular landscape.

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
Whenever there is an image, one can right click and open image in a new tab
to get a higher quality version of the image.

## Files that are ignored

Overall we are ignoring all the rds, rdb and RData files that are generated
in the analysis. We also ignore the cache folders and the files folders. This folder
will be served using a shiny server instead and therefore can be visualized
in a website. The data folders and results folders are also ignored, they would
make the repo too big. Check the gitignore files for a complete list. 
