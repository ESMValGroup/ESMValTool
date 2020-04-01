FROM esmvalgroup/esmvalcore:update_docker

# update the conda packages
RUN conda update -y conda pip

# install development tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    curl \
    unzip

COPY . /source
WORKDIR /source

# install environment packages

RUN conda install -c conda-forge -c esmvalgroup esmvaltool
RUN pip install .
RUN Rscript esmvaltool/install/R/setup.R
RUN julia esmvaltool/install/Julia/setup.jl


# run tests
RUN esmvaltool -h

ENTRYPOINT ["esmvaltool"]
CMD ["-h"]
