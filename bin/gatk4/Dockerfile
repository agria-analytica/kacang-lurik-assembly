FROM mambaorg/micromamba:0.15.3
RUN micromamba install -y -n base gatk4=4.2.3.0 -c bioconda -c conda-forge -c defaults && \
	micromamba clean --all --yes

