FROM mambaorg/micromamba:0.15.3

# install the required compiler to compile swalo.
RUN micromamba install -y -n base -c bioconda -c conda-forge -c defaults \
	compilers \
	c-compiler \
	cxx-compiler \
	gcc_linux-64 \
	gxx_linux-64 \
	make \
	git \
	bowtie2 && \
	micromamba clean --all --yes
RUN git clone https://github.com/atifrahman/SWALO.git && cd SWALO && \
	make && \
	cp {align,bowtie2convert,bowtieconvert,swalo,swaloFile} /opt/conda/bin && \
	cd .. && rm -r SWALO
