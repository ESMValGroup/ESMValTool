FROM continuumio/miniconda3

# update the conda packages
RUN conda update -y conda pip

# install development tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    curl \
    unzip

# create directory for downloaded files
RUN mkdir /src
WORKDIR /src

# install environment packages
RUN conda install -c conda-forge -c esmvalgroup -c birdhouse esmvaltool

# run tests
RUN esmvaltool -h

ENTRYPOINT ["esmvaltool"]
CMD ["-h"]
