# To build this container, go to ESMValTool root folder and execute:
# docker build -t esmvaltool:latest . -f docker/Dockerfile
FROM condaforge/mambaforge

WORKDIR /src/ESMValTool
COPY environment.yml .
RUN mamba update -y conda mamba pip && mamba env create --name esmvaltool --file environment.yml && conda clean --all -y

# Make RUN commands use the new environment:
SHELL ["conda", "run", "--name", "esmvaltool", "/bin/bash", "-c"]

COPY . .
RUN pip install --no-cache . && esmvaltool install Julia

ENTRYPOINT ["conda", "run", "--name", "esmvaltool", "esmvaltool"]
