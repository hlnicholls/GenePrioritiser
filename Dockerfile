# Use Miniconda3 base image
FROM continuumio/miniconda3

# Set working directory inside container
WORKDIR /GenePrioritiser

# Copy everything from your local repo into the image
COPY . .

# Create the Conda environment from the exported YAML
RUN conda env create -f GenePrioritiser_env.yml

# Ensure the environment is activated for all RUN/CMD after this
SHELL ["conda", "run", "--no-capture-output", "-n", "GenePrioritiser_env", "/bin/bash", "-c"]

# Optional: install Nextflow if not part of your environment
RUN curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/

# Set default command to run the pipeline
CMD ["nextflow", "run", "main.nf"]

# Install R separately (outside Conda env)
RUN apt-get update && apt-get install -y r-base

# Install required R packages
RUN Rscript -e "install.packages(c('tidyverse', 'magrittr', 'data.table', 'enrichR', 'org.Hs.eg.db', 'tidygraph', 'plyr', 'janitor'), repos='https://cloud.r-project.org')"

