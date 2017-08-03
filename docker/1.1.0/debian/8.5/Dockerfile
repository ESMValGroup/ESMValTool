FROM continuumio/miniconda

# install development tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    curl \
    unzip

# create directory for downloaded files
RUN mkdir /src
WORKDIR /src

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
