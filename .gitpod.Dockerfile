FROM gitpod/workspace-base@sha256:12853f7c901eb2b677a549cb112c85f9679d18feb30093bcc63aa252540ecad9

USER root
ENV TZ=America/Los_Angeles

# Install system dependencies
RUN apt-get update --quiet && \
    apt-get install --quiet --yes --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    wget \
    make \
    texlive-full \
    g++ \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniforge3-Linux-x86_64.sh

# Set PATH for Conda
ENV PATH="/opt/conda/bin:$PATH"

# Change ownership for gitpod
RUN chown -R gitpod:gitpod /opt/conda

# Change user to gitpod
USER gitpod

# Configure conda
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda install --quiet --yes --update-all --name base mamba pandoc && \
    conda clean --all --force-pkgs-dirs --yes
# Install mamba and pandoc
RUN conda install -n base -c conda-forge mamba pandoc -y

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN mamba env create -f /tmp/environment.yml \
    && conda clean -afy

# Enable Jupyter extensions
RUN /bin/bash -c "source activate gl4u_rnaseq_2024 \
    && jupyter contrib nbextension install --user \
    && jupyter nbextensions_configurator enable --user"

# Install RSEM from source
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar -xzvf v1.3.3.tar.gz && \
    cd RSEM-1.3.3 && \
    make && \
    sudo make install && \
    cd .. && \
    rm -rf RSEM-1.3.3 v1.3.3.tar.gz

# Add RSEM to PATH
ENV PATH="/usr/local/bin:${PATH}"

# Configure R
RUN mkdir -p ~/.R && \
    echo "options(repos = c(CRAN = 'https://cloud.r-project.org'))" > ~/.Rprofile

# Install additional R packages
RUN conda run -n gl4u_rnaseq_2024 R -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', dependencies = TRUE); \
    BiocManager::install(c( \
        'tximport', \
        'DESeq2', \
        'org.Mm.eg.db', \
        'org.At.tair.db', \
        'org.Ce.eg.db', \
        'org.Dr.eg.db', \
        'org.Dm.eg.db', \
        'org.Hs.eg.db', \
        'org.Rn.eg.db', \
        'org.Sc.sgd.db', \
        'STRINGdb', \
        'PANTHER.db', \
        'ComplexHeatmap', \
        'EnhancedVolcano', \
        'clusterProfiler', \
        'goseq', \
        'fgsea', \
        'enrichplot' \
    ), ask = FALSE); \
    install.packages('tidyHeatmap', dependencies = TRUE);"

# Install RSeQC 5.0.3 and add to PATH
RUN conda run -n gl4u_rnaseq_2024 pip install RSeQC==5.0.3 && \
    echo 'export PATH=$PATH:$HOME/miniconda3/envs/gl4u_rnaseq_2024/bin' >> $HOME/.bashrc
