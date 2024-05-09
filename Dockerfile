# Use nf-core base image
FROM nfcore/base:latest

LABEL authors="Urmo Võsa" \
      description="Docker image containing all requirements ABF finemapping"

# Copy your environment.yml into the Docker image
COPY environment.yml /
RUN apt-get update && apt install -y libgmp-dev && apt install -y build-essential

# Use Conda to create the environment as per the environment.yml file
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/abfdeps/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
