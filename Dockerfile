FROM mambaorg/micromamba:2.4.0

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /app/

RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

COPY main.nf nextflow.config /app/

# Re-use the base image's entrypoint script
# We must prepend "/usr/local/bin/_entrypoint.sh" to ensure the env is activated
# before nextflow runs.
# ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "nextflow", "run", "main.nf"]