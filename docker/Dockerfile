FROM esmvalgroup/esmvalcore:latest

COPY . /esmvaltool
WORKDIR /esmvaltool

# install environment packages
RUN apt install julia -y
RUN conda env update --name base --file environment.yml && conda clean --all -y
RUN pip install .
RUN Rscript esmvaltool/install/R/setup.R
RUN julia esmvaltool/install/Julia/setup.jl

# run tests
RUN esmvaltool version

ENTRYPOINT ["esmvaltool"]
CMD ["-h"]
