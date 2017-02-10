FROM centos:7

# install development tools
RUN yum clean all && yum makecache && yum update -y && \
    yum groupinstall -y "Development tools"

# create directory for downloaded files
RUN mkdir /src
WORKDIR /src

# isntall conda
RUN curl -L -O https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN bash ./Miniconda2-latest-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda update -y pip

# install python packages specified in conda environment file
COPY environment.yml /src/environment.yml
RUN conda env update -f environment.yml

# install esmvaltool
RUN curl -L -O https://github.com/ESMValGroup/ESMValTool/archive/docker.zip
RUN unzip docker.zip
WORKDIR /src/ESMValTool-docker
COPY config_private.xml /src/ESMValTool-docker/config_private.xml

ENTRYPOINT ["python", "main.py"]
CMD ["nml/namelist_test.xml"]
