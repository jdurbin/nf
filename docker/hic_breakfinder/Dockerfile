FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive


# hic_breakfinder installation is pretty fragile. 
# It needs JasonCpp, Eigen, and bamtools installed first.  
# Note that bamtools must be version 2.4.2.   Not 2.4.0, or 2.4.1, or 2.5.x. 
# It's also very picky about setting CPPFLAGS
# Couldn't get configure to add -lz properly, so ended up just doing the linking step manually (g++ ...)
# This also means the make will fail generating an error, thus "make; exit 0"
# None of this is ideal and it'd be good to figure out how to get ./configure;make to work properly
# but this at least gives us something that works and documents the issues for whoever wants to take that on. 

RUN apt-get update -y\    
	&& DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y \
    awscli \
	python3 \
	python3-pip \
	git \
	libbz2-dev \
	lz4 \
	zlib1g-dev \
	build-essential \
	autotools-dev \
	automake \
    wget \
    cmake \
    libjsoncpp-dev \
    && wget https://gitlab.com/libeigen/eigen/-/archive/3.4-rc1/eigen-3.4-rc1.tar.gz \
    && tar -zxvf eigen-3.4-rc1.tar.gz \
    && cp -r eigen-3.4-rc1/Eigen /usr/local/include/ \
    && cd / \
    && wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.4.2.tar.gz \
    && tar -zxvf v2.4.2.tar.gz \
    && cd bamtools-2.4.2 \ 
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \ 
    && cd / \       
    && git clone https://github.com/dixonlab/hic_breakfinder.git \
    && cd hic_breakfinder \
    && ./configure CPPFLAGS="-I /usr/local/include/bamtools -I /usr/local/include/Eigen/" LDFLAGS="-L /usr/local/lib -lz" \
    && make; exit 0 \
    && cd src \
    && g++  -g -O2 -g -L /usr/local/lib -L /lib/x86_64-linux-gnu/ -o hic_breakfinder hic_breakfinder.o  -lbamtools -lz \
    && cd .. \
    && make install \	
	&& apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*






