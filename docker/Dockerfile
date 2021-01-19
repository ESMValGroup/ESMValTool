# To build this container, go to ESMValTool root folder and execute:
# docker build -t esmvaltool:latest . -f docker/Dockerfile
FROM continuumio/miniconda3

RUN apt install julia -y && apt-get clean
COPY ./environment.yml /src/ESMValTool/environment.yml
WORKDIR /src/ESMValTool
RUN conda update -y conda && conda env create --name esmvaltool --file environment.yml && conda clean --all -y

# Make RUN commands use the new environment:
SHELL ["conda", "run", "--name", "esmvaltool", "/bin/bash", "-c"]

COPY ./esmvaltool/install/R /src/ESMValTool/esmvaltool/install/R
RUN Rscript ./esmvaltool/install/R/setup.R

COPY . /src/ESMValTool
RUN pip install -e . && pip cache purge && esmvaltool install Julia

ENTRYPOINT ["conda", "run", "--name", "esmvaltool", "esmvaltool"]
