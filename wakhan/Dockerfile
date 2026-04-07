FROM condaforge/mambaforge:latest

LABEL maintainer="Wakhan Docker"

# Install mamba for faster package resolution
RUN mamba update -n base -c defaults mamba -y

# Create conda environment with wakhan
RUN mamba create -n wakhan_env -c conda-forge -c bioconda -c defaults \
    wakhan=0.4.2 \
    -y && \
    mamba clean --all -y

# Make wakhan available in PATH
ENV PATH="/opt/conda/envs/wakhan_env/bin:${PATH}"

# Activate conda environment by default
SHELL ["/bin/bash", "-c"]
RUN echo "source activate wakhan_env" > ~/.bashrc

WORKDIR /data

CMD ["/bin/bash"]
