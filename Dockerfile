# ============================
# Stage 1 — Build
# ============================
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

# Install build dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
        wget \
        git \
        autoconf \
        automake \
        libtool \
        autogen \
        build-essential \
        curl \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libnss-sss \
        libssl-dev \
        libxml2-dev \
        libuv1-dev \
        zlib1g-dev \
        ncbi-blast+ \
        r-base \
        r-bioc-rsamtools \
        r-bioc-biostrings \
        r-cran-optparse \
        r-cran-stringr && \
    rm -rf /var/lib/apt/lists/*

# Install rBLAST via devtools, same as SCRAMble
RUN Rscript -e " \
    install.packages('devtools', repos='https://cloud.r-project.org'); \
    install.packages('remotes', repos='https://cloud.r-project.org'); \
    remotes::install_github('mhahsler/rBLAST') \
"

# Build htslib from source
RUN git clone --branch 1.21 --recurse-submodules https://github.com/samtools/htslib.git && \
    cd htslib && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install && \
    echo "/usr/local/lib" > /etc/ld.so.conf.d/local.conf && ldconfig && \
    cd .. && rm -rf htslib

# Build bcftools from source (same version as htslib)
RUN git clone --branch 1.21 https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    autoheader && autoconf && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf bcftools

# Install Miniforge — uses conda-forge by default, no Anaconda ToS issue
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

# Install python dependencies
RUN /opt/conda/bin/conda install -y \
        python=3.11 \
        pandas && \
    /opt/conda/bin/conda clean -afy && \
    pip install pysam dxpy


# Copy source code and build scramble
RUN mkdir -p /app
COPY . /app
RUN cd /app/cluster_identifier/src && \
    make clean && \
    make


# ============================
# Stage 2 — Runtime
# ============================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

# Install runtime system dependencies only
RUN apt-get update && apt-get install -y \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libncurses5 \
        libssl3 \
        libxml2 \
        zlib1g \
        libnss-sss \
        ncbi-blast+ \
        libuv1-dev \
        r-base \
        r-bioc-rsamtools \
        r-bioc-biostrings \
        r-cran-optparse \
        r-cran-stringr && \
    rm -rf /var/lib/apt/lists/*

# Copy htslib shared libs from builder
COPY --from=builder /usr/local/lib/libhts.so* /usr/local/lib/
# Copy bcftools and bgzip binaries from builder
COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/
COPY --from=builder /usr/local/bin/bgzip /usr/local/bin/
RUN echo "/usr/local/lib" > /etc/ld.so.conf.d/local.conf && ldconfig

# Copy Miniforge from builder — includes R, all R packages, Python, pysam, pandas
COPY --from=builder /opt/conda /opt/conda

# Copy R packages installed in builder (rBLAST)
COPY --from=builder /usr/local/lib/R /usr/local/lib/R

# Copy compiled cluster_identifier binary
COPY --from=builder /app/cluster_identifier/src/build/cluster_identifier /usr/local/bin/

# Copy app folder (scramble R scripts + Python scripts)
COPY --from=builder /app /app

# Create output directory
RUN mkdir -p /app/output

ENTRYPOINT ["python", "/app/ALU_analysis.py"]
