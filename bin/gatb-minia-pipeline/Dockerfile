FROM mambaorg/micromamba:0.15.3
RUN micromamba install -y -n base -c bioconda -c conda-forge -c defaults \
	git \
	python=2.7 \
	scipy \
	numpy \
	mathstats \
	pysam=0.8.3 \
	bwa \
	samtools=1.13 \
	networkx=1.11 \
	matplotlib && \
	micromamba clean --all --yes
RUN git clone --recursive https://github.com/GATB/gatb-minia-pipeline .

