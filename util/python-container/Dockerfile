FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies for scientific/bioinformatics Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libffi-dev \
    gfortran \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    # Data science core
    numpy \
    pandas \
    scipy \
    scikit-learn \
    statsmodels \
    # Visualization
    matplotlib \
    seaborn \
    plotly \
    # Bioinformatics
    pysam \
    biopython \
    pybedtools \
    pyBigWig \
    # Utilities
    click \
    tqdm \
    pyyaml \
    requests \
    polars \
    pyarrow \
    # File formats
    openpyxl \
    xlsxwriter

# Strip build dependencies to reduce image size
RUN apt-get update && apt-get remove -y \
    build-essential \
    gfortran \
    pkg-config \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Keep runtime libs needed by compiled packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Verify key packages
RUN python3 -c "import numpy; print(f'numpy {numpy.__version__}')" && \
    python3 -c "import pandas; print(f'pandas {pandas.__version__}')" && \
    python3 -c "import pysam; print(f'pysam {pysam.__version__}')" && \
    python3 -c "import matplotlib; print(f'matplotlib {matplotlib.__version__}')"

WORKDIR /data

ENTRYPOINT ["python3"]
