# To build this container, go to ESMValCore root folder and execute :
# docker build -t ${TAG_OF_YOUR_CHOICE} . -f docker/Dockerfile
FROM continuumio/miniconda3

COPY . /esmvaltool
WORKDIR /esmvaltool

RUN conda update -y conda pip && conda env update --name base --file environment.yml && conda clean --all -y
RUN pip install . && pip cache purge
RUN apt install julia -y
RUN julia esmvaltool/install/Julia/setup.jl
RUN Rscript esmvaltool/install/R/setup.R
RUN rm -r /esmvaltool

RUN esmvaltool -v

ENTRYPOINT ["esmvaltool"]
CMD ["-h"]
