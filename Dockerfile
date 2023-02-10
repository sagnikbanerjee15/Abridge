###################################################################################################################################################################################################
# ABRIDGE
###################################################################################################################################################################################################
FROM rust
LABEL maintainer="sagnikbanerjee15@gmail.com" 

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

# Update base image and install software
RUN apt-get -y update
RUN apt-get -y install --fix-missing git python3 less vim wget time zlib1g zlib1g-dev lzma-dev \
	libncurses5-dev libcurl4-nss-dev liblzma-dev libncursesw5-dev make unzip zip build-essential \
	gcc g++ cmake ca-certificates libbz2-dev xz-utils htop autoconf automake binutils bison flex \
	gettext libtool patch pkg-config p7zip-full p7zip pip perl libc6-dbg gdb valgrind
RUN apt-get clean all

# Install base utilities
RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda
     
# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda update -y conda

#######################################################################################################################################################################
# FCLQC
#######################################################################################################################################################################

RUN mkdir -p /software && cd /software
RUN cd /software && git clone https://github.com/Minhyeok01/FCLQC.git
RUN cd /software/FCLQC && cargo build --release

ENV PATH=$PATH:/software/FCLQC/target/release

#######################################################################################################################################################################
# ZPAQ
#######################################################################################################################################################################
RUN conda install -c sbu-hpc -y zpaq=7.15

#######################################################################################################################################################################
# SAMTOOLS
#######################################################################################################################################################################

RUN apt-get update -y && apt-get install -y \
    apt-utils \
    bzip2 \
    gcc \
    make \
    ncurses-dev \
    wget \
    zlib1g-dev

##############
#HTSlib 1.16#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.16.tar.bz2

WORKDIR /tmp/htslib-1.16
RUN ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

################
#Samtools 1.16.1#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar --bzip2 -xf samtools-1.16.1.tar.bz2

WORKDIR /tmp/samtools-1.16.1
RUN ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install

WORKDIR /
RUN rm -rf /tmp/samtools-1.16.1

#######
#tabix#
#######
RUN ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

###################################################################################################################################################################################################
# ABRIDGE
###################################################################################################################################################################################################
RUN pip install ruffus

ARG ABRIDGE_VERSION=1.2.0
COPY marker /dev/null
RUN mkdir -p /software/abridge && cd /software/abridge && \
	git clone --recurse-submodules https://github.com/sagnikbanerjee15/Abridge.git && \
	cd Abridge && \
	make htslib
	

RUN cd /software/abridge/Abridge && \
	git pull && \
	make install
	
ENV PATH=${PATH}:/software/abridge/Abridge:/software/abridge/Abridge/src:/software/abridge/Abridge/bin
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"

ENV PATH=${PATH}:/opt/samtools/bin
