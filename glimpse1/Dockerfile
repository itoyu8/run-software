FROM ubuntu:20.04

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    autoconf \
    autotools-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Install HTSlib v1.11
RUN wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
    tar -xf htslib-1.11.tar.bz2 && \
    mv htslib-1.11 htslib && \
    cd htslib && \
    autoheader && autoconf && ./configure && make && \
    cd .. && rm htslib-1.11.tar.bz2

# Install Boost v1.73.0
RUN wget https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2 && \
    tar --bzip2 -xf boost_1_73_0.tar.bz2 && \
    cd boost_1_73_0 && \
    ./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
    ./b2 install && \
    cd .. && rm boost_1_73_0.tar.bz2

# Set environment variables for compilation
ENV HTSSRC=/opt/htslib
ENV HTSLIB_INC=/opt/htslib
ENV HTSLIB_LIB=/opt/htslib
ENV BOOST_INC=/opt/boost/include
ENV BOOST_LIB_IO=/opt/boost/lib/libboost_iostreams.a
ENV BOOST_LIB_PO=/opt/boost/lib/libboost_program_options.a

# Clone and compile GLIMPSE1 (not GLIMPSE2)
RUN git clone https://github.com/odelaneau/GLIMPSE.git && \
    cd GLIMPSE && \
    git checkout v1.1.1

# Compile all GLIMPSE1 tools  
WORKDIR /opt/GLIMPSE
RUN ls -la && \
    for tool in chunk ligate phase split_reference sample concordance; do \
        echo "Processing GLIMPSE1 tool: $tool"; \
        if [ -d "$tool" ]; then \
            cd $tool && \
            echo "Updating makefile for $tool..." && \
            sed -i 's|HTSSRC=.*|HTSSRC=/opt/htslib|g' makefile && \
            sed -i 's|HTSLIB_INC=.*|HTSLIB_INC=/opt/htslib|g' makefile && \
            sed -i 's|HTSLIB_LIB=.*|HTSLIB_LIB=/opt/htslib/libhts.a|g' makefile && \
            sed -i 's|BOOST_INC=.*|BOOST_INC=/opt/boost/include|g' makefile && \
            sed -i 's|BOOST_LIB_IO=.*|BOOST_LIB_IO=/opt/boost/lib/libboost_iostreams.a|g' makefile && \
            sed -i 's|BOOST_LIB_PO=.*|BOOST_LIB_PO=/opt/boost/lib/libboost_program_options.a|g' makefile && \
            sed -i 's|/usr/lib/x86_64-linux-gnu/libboost_serialization.a|/opt/boost/lib/libboost_serialization.a|g' makefile && \
            make && \
            cd ..; \
        else \
            echo "Directory $tool not found"; \
        fi; \
    done

# Add GLIMPSE1 binaries to PATH
ENV PATH="/opt/GLIMPSE/chunk/bin:/opt/GLIMPSE/ligate/bin:/opt/GLIMPSE/phase/bin:/opt/GLIMPSE/split_reference/bin:/opt/GLIMPSE/sample/bin:/opt/GLIMPSE/concordance/bin:${PATH}"

# Set final working directory
WORKDIR /data

CMD ["/bin/bash"]