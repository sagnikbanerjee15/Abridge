###################################################################################################################################################################################################
# ABRIDGE
###################################################################################################################################################################################################
FROM ubuntu:22.04
LABEL maintainer="sagnikbanerjee15@gmail.com" 
LABEL org.opencontainers.image.source https://github.com/sagnikbanerjee15/dockerized_tools_and_pipelines

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

# Update base image and install software
RUN apt-get -y update
RUN apt-get -y install --fix-missing git python3 less vim wget time zlib1g zlib1g-dev lzma-dev \
	libncurses5-dev libcurl4-nss-dev liblzma-dev libncursesw5-dev make unzip zip build-essential \
	gcc g++ cmake ca-certificates libbz2-dev xz-utils htop autoconf automake binutils bison flex \
	gettext libtool patch pkg-config p7zip-full p7zip pip perl libc6-dbg gdb valgrind
RUN apt-get clean all

###################################################################################################################################################################################################
# ABRIDGE
###################################################################################################################################################################################################
RUN pip install ruffus
COPY marker /dev/null
ARG ABRIDGE_VERSION=1.1.0
RUN mkdir -p /software/abridge && cd /software/abridge && \
	git clone --recurse-submodules https://github.com/sagnikbanerjee15/Abridge.git && \
	cd Abridge && \
	make htslib && \
	make install
	
ENV PATH=${PATH}:/software/abridge/Abridge:/software/abridge/Abridge/src:/software/abridge/Abridge/bin:~/work/ABRIDGE/Abridge/src:~/work/ABRIDGE/Abridge:~/work/ABRIDGE/Abridge/bin