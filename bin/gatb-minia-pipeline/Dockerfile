FROM mambaorg/micromamba:0.15.3
RUN micromamba install -y -n base -c bioconda -c conda-forge -c default \
	git \
	python=2.7 \
	scipy \
	numpy \
	mathstats \
	pysam \
	bwa && \
	micromamba clean --all --yes
RUN git clone --recursive https://github.com/GATB/gatb-minia-pipeline .

