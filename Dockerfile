# name of the docker image: molecular_landscape
# date: 2023-01-18 (YYYY-MM-DD) 
# run the command to build the docker:
# sudo docker build --cpuset-cpus=0-30 -t molecular_landscape .
# The number of cpus used to get the docker image is 30 

FROM rocker/rstudio:4.2.0

# system libraries of general use
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
        libz-dev \
        zlib1g \
        zlib1g-dev \
        libpng-dev \
        pkg-config \
        libfontconfig1-dev \
        cmake \
        libnlopt-dev \
        libhdf5-serial-dev \
        libhdf5-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libxml2-dev \
        libpng-dev \
        libxt-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libglpk40 \
        libgit2-28 \
        vim \
    && apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install quarto so we can compile the book
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.0.37/quarto-1.0.37-linux-amd64.deb
RUN dpkg -i quarto-1.0.37-linux-amd64.deb
RUN rm quarto-1.0.37-linux-amd64.deb

WORKDIR /home/rstudio/molecular_landscape

# Move the packages file so we can use it to indicate which packages
# will be installed
COPY --chown=rstudio:rstudio renv.lock .Rprofile .gitignore ./
#COPY renv.lock .Rprofile .gitignore ./
COPY --chown=rstudio:rstudio renv renv
#COPY renv renv

# this will be removed and fetched automatically from github
COPY --chown=rstudio:rstudio scripts scripts
#COPY scripts scripts

# we add the R project file because we can then use renv directly
COPY --chown=rstudio:rstudio 20220721_molecular_landscape.Rproj ./
#COPY 20220721_molecular_landscape.Rproj ./
COPY --chown=rstudio:rstudio data data
#COPY data data

RUN bash -c 'mkdir -p \
    results/{plots,rds_files,tables}/{pca_merging,scoring,surv_analysis_estrogen,validation}'
RUN mkdir -p results/book
RUN chown -R rstudio results/

# gives sudo access to everything in this folder, so the user can actually
# create new folders
RUN chown -R rstudio /home/rstudio/molecular_landscape

# initialize the renv repo
RUN R -q -e 'renv::restore()'