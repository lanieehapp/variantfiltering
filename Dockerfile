# Base Image
FROM r-base:4.0.2

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update && apt-get install -y \
    build-essential \
    #libcurl4-gnutls-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \ 
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    gawk \
    unzip \
    less \
    git \
    bedtools \
    parallel \
    samtools 	
# install Python libraries
WORKDIR /usr/local/bin


# install R required dependencies
RUN R --vanilla -e 'install.packages(c("devtools", "stringr", "plyr","doParallel", "data.table", "bedr"), repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'devtools::install_github(repo="knausb/vcfR")'


# clone variantfiltering repo
ADD https://api.github.com/repos/lanieehapp/variantfiltering/git/refs/heads/ version.json
RUN git clone https://github.com/lanieehapp/variantfiltering.git

# add variantfiltering repo to SYSPATH
ENV PATH variantfiltering:$PATH

# change the permission of the repo
RUN chmod 777 -R variantfiltering
WORKDIR /usr/local/bin/variantfiltering

