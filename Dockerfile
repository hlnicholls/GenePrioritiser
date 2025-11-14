# Use Miniconda3 base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /GenePrioritiser

# Copy local code
COPY . .

# Install updated libstdc++ and Java and R
RUN apt-get update && \
    apt-get install -y openjdk-17-jdk r-base build-essential software-properties-common && \
    apt-get install -y libstdc++6 && \
    strings /usr/lib/aarch64-linux-gnu/libstdc++.so.6 | grep GLIBCXX

# Create the Conda environment
RUN conda env create -f GenePrioritiser_env.yml

# Pip installs inside the environment
RUN conda run -n GenePrioritiser_env pip install --force-reinstall scikit-learn==1.4.2 && \
    conda run -n GenePrioritiser_env pip install --force-reinstall scipy==1.11.4 && \
    conda run -n GenePrioritiser_env pip install --force-reinstall numpy==1.23.0 && \
    conda run -n GenePrioritiser_env pip install scikit-optimize

# Install required R packages
#RUN Rscript -e "install.packages(c('tidyverse', 'magrittr', 'data.table', 'enrichR', 'org.Hs.eg.db', 'tidygraph', 'plyr', 'janitor'), repos='https://cloud.r-project.org')"

# Set default command
CMD ["nextflow", "run", "main.nf"]
