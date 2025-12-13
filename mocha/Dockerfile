FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    git \
    g++ \
    make \
    zlib1g-dev \
    samtools \
    bedtools \
    libbz2-dev \
    libssl-dev \
    liblzma-dev \
    libgsl0-dev \
    libcurl4-openssl-dev \
    openjdk-11-jre-headless \
    r-base \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /tmp

# Download and compile BCFtools 1.20 with HTSlib
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.20 bcftools-1.20.tar.bz2

# Download MoChA plugins and compile
RUN cd /tmp && \
    wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    rm -f plugins/{mocha,beta_binom,genome_rules}.h plugins/{mocha,mochatools,extendFMT}.c && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/mocha.h && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/beta_binom.h && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/genome_rules.h && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/mocha.c && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/mochatools.c && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/extendFMT.c && \
    make && \
    mkdir -p /usr/local/libexec/bcftools && \
    cp plugins/fill-tags.so plugins/fixploidy.so plugins/mocha.so plugins/mochatools.so plugins/extendFMT.so /usr/local/libexec/bcftools/ && \
    cd .. && \
    rm -rf bcftools-1.20 bcftools-1.20.tar.bz2

# Set environment variables for BCFtools plugins
ENV BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

# Download and install IMPUTE5
RUN wget -O impute5_v1.2.0.zip https://www.dropbox.com/sh/mwnceyhir8yze2j/AABKBCgZsQqz8TlZGo7yXwx6a/impute5_v1.2.0.zip?dl=0 && \
    unzip -ojd /usr/local/bin impute5_v1.2.0.zip impute5_v1.2.0/impute5_v1.2.0_static impute5_v1.2.0/xcftools_static && \
    chmod a+x /usr/local/bin/impute5_v1.2.0_static /usr/local/bin/xcftools_static && \
    ln -s /usr/local/bin/impute5_v1.2.0_static /usr/local/bin/impute5 && \
    rm impute5_v1.2.0.zip

# Download and install Beagle
RUN wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar -O /usr/local/bin/beagle.jar

# Download MoChA R scripts
RUN mkdir -p /usr/local/share/mocha && \
    wget -P /usr/local/share/mocha https://raw.githubusercontent.com/freeseek/mocha/master/mocha_plot.R && \
    wget -P /usr/local/share/mocha https://raw.githubusercontent.com/freeseek/mocha/master/summary_plot.R && \
    wget -P /usr/local/share/mocha https://raw.githubusercontent.com/freeseek/mocha/master/pileup_plot.R && \
    chmod a+x /usr/local/share/mocha/*.R

WORKDIR /data

CMD ["/bin/bash"]
