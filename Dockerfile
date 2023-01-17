# name of the docker image: molecular_landscape
# date: 2023-01-13 (YYYY-MM-DD) 
# run the command to build the docker:
# sudo docker build --cpuset-cpus=0-30 -t molecular_landscape .
# The number of cpus used to get the docker image is 30 

FROM rocker/rstudio:4.2.0

WORKDIR /home/rstudio/molecular_landscape

# Move the packages file so we can use it to indicate which packages
# will be installed
COPY renv.lock .Rprofile .gitignore ./
COPY renv renv

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
    && apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
    
RUN R -q -e 'renv::restore()'

COPY data data
