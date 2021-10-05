FROM mambaorg/micromamba:0.15.3
RUN micromamba install -y -n base -c bioconda -c conda-forge -c default \
	quast=5.0.2 && \
	micromamba clean --all --yes

