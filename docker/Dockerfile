# Recipe from: https://tech.quantco.com/blog/pixi-production
# To build this container, go to ESMValTool root folder and execute:
# docker build --pull -t esmvaltool:latest . -f docker/Dockerfile
FROM ghcr.io/prefix-dev/pixi:latest AS build

WORKDIR /app
COPY . .
RUN pixi global install git
RUN pixi install --frozen -e esmvaltool
RUN pixi shell-hook -e esmvaltool -s bash > /shell-hook
RUN echo "#!/bin/bash" > /app/entrypoint.sh
RUN cat /shell-hook >> /app/entrypoint.sh
RUN echo 'esmvaltool "$@"' >> /app/entrypoint.sh

FROM ubuntu:latest AS production
WORKDIR /app
# Copy the environment:
COPY --from=build /app/.pixi/envs/esmvaltool /app/.pixi/envs/esmvaltool
# Copy the source code because esmvaltool is installed in editable mode:
COPY --from=build /app/esmvaltool /app/esmvaltool
# Copy the entrypoint script and make it executable:
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh

ENTRYPOINT [ "/app/entrypoint.sh" ]
