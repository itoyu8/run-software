FROM condaforge/mambaforge:latest

LABEL maintainer="Severus Docker"

# Install mamba for faster package resolution
RUN mamba update -n base -c defaults mamba -y

# Create conda environment with severus
RUN mamba create -n severus_env -c conda-forge -c bioconda -c defaults \
    python=3.10 \
    samtools>=1.14 \
    networkx>=2.6 \
    pygraphviz \
    pydot \
    matplotlib-base \
    biopython \
    numpy \
    pysam \
    plotly \
    -y && \
    mamba clean --all -y

# Install severus from GitHub
RUN /bin/bash -c "source activate severus_env && \
    pip install git+https://github.com/KolmogorovLab/Severus.git"

# Make severus available in PATH
ENV PATH="/opt/conda/envs/severus_env/bin:${PATH}"

# Activate conda environment by default
SHELL ["/bin/bash", "-c"]
RUN echo "source activate severus_env" > ~/.bashrc

WORKDIR /data

CMD ["/bin/bash"]
