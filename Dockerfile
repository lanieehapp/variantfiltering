# Base Image
FROM r-base:3.5.2

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    curl \
    vim \
    devscripts \
    less \
    wget \
    unzip \
    cmake \
    python3 \
    gawk \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    autoconf \
    bsdmainutils \
    bedtools \
    gcc-8-base \
    libmpx2 \		
    libgcc-8-dev \
    bedops \
    tabix \
    parallel 	
# install Python libraries
WORKDIR /usr/local/bin


# install R required dependencies
RUN R --vanilla -e 'install.packages(c("devtools", "stringr", "plyr","doParallel", "data.table", "qqman", "bedr"), repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'devtools::install_github(repo="knausb/vcfR")'


# clone disco-wave repo
ADD https://api.github.com/repos/lanieehapp/variantfiltering/git/refs/heads/ version.json
RUN git clone https://github.com/lanieehapp/variantfiltering.git

# add disco-wave repo to SYSPATH
ENV PATH variantfiltering:$PATH

# change the permission of the repo
RUN chmod 777 -R variantfiltering
WORKDIR /usr/local/bin/variantfiltering
