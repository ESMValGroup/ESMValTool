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

# install python packages specified in conda environment file (copied from context)
COPY environment.yml /src/environment.yml
RUN conda env update -f /src/environment.yml

#Copy entire ESMValTool source into the container
RUN mkdir /src/ESMValTool
COPY . /src/ESMValTool/

WORKDIR /src/ESMValTool

#overwrite default config_private with one specifically created for Docker
COPY docker/1.1.0/centos/7/config_private.xml /src/ESMValTool/config_private.xml

ENTRYPOINT ["python", "main.py"]
CMD ["nml/namelist_test.xml"]
