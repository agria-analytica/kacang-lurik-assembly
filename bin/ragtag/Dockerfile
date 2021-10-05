FROM mambaorg/micromamba:0.15.3
RUN micromamba install -y -n base -c bioconda -c conda-forge -c defaults \
	ragtag=2.0.1 && \
	micromamba clean --all --yes

