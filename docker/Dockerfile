# To build this container, go to ESMValTool root folder and execute:
# docker build -t esmvaltool:latest . -f docker/Dockerfile
FROM continuumio/miniconda3

COPY ./conda-linux-64.lock /src/ESMValTool/conda-linux-64.lock
WORKDIR /src/ESMValTool
RUN conda install -y -c conda-forge mamba && mamba create --name esmvaltool --file conda-linux-64.lock && conda clean --all -y

# Make RUN commands use the new environment:
SHELL ["conda", "run", "--name", "esmvaltool", "/bin/bash", "-c"]

COPY ./esmvaltool/install/R /src/ESMValTool/esmvaltool/install/R
RUN Rscript ./esmvaltool/install/R/setup.R

COPY . /src/ESMValTool
RUN pip install -e . && pip cache purge && esmvaltool install Julia

ENTRYPOINT ["conda", "run", "--name", "esmvaltool", "esmvaltool"]
