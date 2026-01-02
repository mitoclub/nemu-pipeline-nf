# Use the mamba image as the base image
FROM mambaorg/micromamba:2.4.0

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY environment.yml /app/environment.yml

# Create the conda environment using mamba  TODO --use-uv
RUN micromamba create -f environment.yml --name nemu-pipeline --yes && \
    micromamba clean --all --yes

# Activate the environment
SHELL ["/bin/bash", "-c"]
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate nemu-pipeline" >> ~/.bashrc

# Copy the pipeline scripts into the container
COPY main.nf /app

# Set the entry point to the Nextflow script
ENTRYPOINT ["nextflow", "run", "main.nf"]