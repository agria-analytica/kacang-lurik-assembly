# Genome assembly Arachis hypogea var. kacang lurik.
SHELL = /bin/bash
DOCKER = docker run --rm -v $$(pwd):/project -w /project
DOCKER_FASTX = $(DOCKER) quay.io/biocontainers/fastx_toolkit:0.0.14--he1b5a44_8
DOCKER_FASTQC = $(DOCKER) quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1
DOCKER_KMC = $(DOCKER) quay.io/biocontainers/kmc:3.1.2rc1--h2d02072_0
DOCKER_TRIMMOMATIC = $(DOCKER) quay.io/gitobioinformatics/trimmomatic:0.38

DIRS = results results/fastx results/fastqc \
       results/kmc results/kmc/tmp \
	   results/trimmomatic
DATA = data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_1.fq.gz \
       data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_2.fq.gz

all: $(DIRS) qc trimmomatic genomescope

$(DIRS): 
	[ -d $@ ] || mkdir $@

qc: fastx fastqc

FASTX_RESULTS = $(addprefix results/fastx/, $(addsuffix .stats.txt, $(basename $(notdir $(DATA)))))
fastx: $(FASTX_RESULTS)

$(FASTX_RESULTS): results/fastx/%.stats.txt: data/X401SC21062291-Z01-F001/raw_data/KL/%.gz
	$(DOCKER_FASTX) sh -c "zcat $< | fastx_quality_stats -o $@" && \
	$(DOCKER_FASTX) fastq_quality_boxplot_graph.sh -i $@ -o $(basename $@).qual.png && \
	$(DOCKER_FASTX) fastx_nucleotide_distribution_graph.sh -i $@ -o $(basename $@).nuc.png

FASTQC_RESULTS = $(addprefix results/fastqc/, $(addsuffix _fastqc.html, $(basename $(basename $(notdir $(DATA))))))
fastqc: $(FASTQC_RESULTS)

$(FASTQC_RESULTS): results/fastqc/%_fastqc.html: data/X401SC21062291-Z01-F001/raw_data/KL/%.fq.gz
	echo $@; \
	$(DOCKER_FASTQC) fastqc -t 4 -o results/fastqc  $<

# trimmomatic.
TRIM_FASTQ_PAIRED1 = $(addprefix results/trimmomatic/, $(addsuffix .trim.gz, $(basename $(notdir data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_1.fq.gz))))
TRIM_FASTQ_PAIRED2 = $(addprefix results/trimmomatic/, $(addsuffix .trim.gz, $(basename $(notdir data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_2.fq.gz))))
TRIM_FASTQ_UNPAIRED1 = $(addsuffix .unpaired.gz, $(basename $(TRIM_FASTQ_PAIRED1)))
TRIM_FASTQ_UNPAIRED2 =  $(addsuffix .unpaired.gz, $(basename $(TRIM_FASTQ_PAIRED2)))
ADAPTER = illumina-adapter-novogene.fa
LOG_TRIMMOMATIC = results/trimmomatic/log.out
trimmomatic: $(TRIM_FASTQ_PAIRED1)

$(TRIM_FASTQ_PAIRED1): $(DATA)
	$(DOCKER_TRIMMOMATIC) PE -threads 10 -trimlog $(LOG_TRIMMOMATIC) $^ $@ $(TRIM_FASTQ_UNPAIRED1) $(TRIM_FASTQ_PAIRED2) \
	$(TRIM_FASTQ_UNPAIRED2) ILLUMINACLIP:$(ADAPTER):2:30:10 SLIDINGWINDOW:3:5

# genomescope
KMER_TABLE = $(addsuffix .kmc_pre, $(addprefix results/kmc/, $(notdir $(basename $(DATA)))))
KMER_HIST  = results/kmc/kacang.lurik.hist
KMER_MERGE  = results/kmc/kacang.lurik.merge.kmc_pre
genomescope: kmc

kmc: $(KMER_HIST)

$(KMER_HIST): $(KMER_MERGE)
	$(DOCKER_KMC) kmc_tools transform $< histogram $@ -cx10000

# merge kmer from all fastq.
$(KMER_MERGE): $(KMER_TABLE)
	$(DOCKER_KMC) kmc_tools simple $(basename $^) union $(basename $@)

# generate kmer tabel from raw data.
$(KMER_TABLE): results/kmc/%.kmc_pre: data/X401SC21062291-Z01-F001/raw_data/KL/%.gz
	$(DOCKER_KMC) kmc -k21 -t10 -m48 -ci3 -cs10000 $< $(basename $@) results/kmc/tmp/


