FROM nvidia/cuda:8.0-cudnn5-devel-ubuntu16.04
LABEL maintainer dlindenbaum


## Install General Requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
        apt-utils \
        build-essential \
        cmake \
        git \
        wget \
        vim \
        python-dev \
        python-pip \
        python-setuptools
RUN pip install --upgrade pip

## Install Basics for Python
RUN apt-get update && apt-get install -y --no-install-recommends \
        python-numpy \
        python-scipy

## Install Essentials for Pillow
RUN apt-get update && apt-get install -y --no-install-recommends \
        libjpeg-dev \
        zlib1g \
        libtiff5-dev \
        libfreetype6-dev \
        libwebp-dev \
        libopenjpeg-dev

RUN pip install Pillow

## Install GDAL Requirments
RUN apt-get update && apt-get install -y --no-install-recommends \
        gdal-bin \
        python-gdal

## Instal OpenCV Requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
        libopencv-dev \
        python-opencv

## Install RTRee
RUN apt-get update && apt-get install -y --no-install-recommends \
        libspatialindex-dev

RUN pip install rtree

ENV GIT_BASE=/opt/
WORKDIR $GIT_BASE

# Download spaceNetUtilities
# FIXME: use ARG instead of ENV once DockerHub supports this

RUN git clone --depth 1 https://github.com/SpaceNetChallenge/utilities.git

ENV PYUTILS_ROOT $GIT_BASE/utilities/python
ENV PYTHONPATH $PYUTILS_ROOT:$PYTHONPATH

WORKDIR /workspace


