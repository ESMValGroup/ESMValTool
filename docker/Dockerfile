# To build this container, go to ESMValCore root folder and execute :
# docker build -t ${TAG_OF_YOUR_CHOICE} . -f docker/Dockerfile
FROM continuumio/miniconda3

RUN apt install julia tree -y && apt-get clean
COPY ./environment.yml /src/ESMValTool/environment.yml
WORKDIR /src/ESMValTool
RUN conda update -y conda pip && conda env update --name base --file environment.yml && conda clean --all -y

COPY ./esmvaltool/install/R /src/ESMValTool/esmvaltool/install/R
RUN Rscript ./esmvaltool/install/R/setup.R

COPY . /src/ESMValTool
RUN pip install . && pip cache purge && esmvaltool --help && esmvaltool install Julia
RUN rm -rf /src


ENTRYPOINT ["esmvaltool"]
