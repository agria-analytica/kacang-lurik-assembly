# Genome assembly Arachis hypogea var. kacang lurik.
SHELL = /bin/bash
DOCKER = docker run --rm -v $$(pwd):/project -w /project
DOCKER_FASTX = $(DOCKER) quay.io/biocontainers/fastx_toolkit:0.0.14--he1b5a44_8
DOCKER_FASTQC = $(DOCKER) quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1
DOCKER_KMC = $(DOCKER) quay.io/biocontainers/kmc:3.1.2rc1--h2d02072_0
DOCKER_TRIMMOMATIC = $(DOCKER) quay.io/gitobioinformatics/trimmomatic:0.38
DOCKER_GATB_PIPELINE = docker run --rm -v $$(pwd):/tmp/project -w /tmp/project agria.analytica/gatb-minia-pipeline:9d56f42
DOCKER_QUAST = $(DOCKER) -v $$(pwd)/../ref-genomes/Arhypogaea/var.Tifrunner:/project/genome agria.analytica/quast:5.0.2
DOCKER_RAGTAG = $(DOCKER) agria.analytica/ragtag:2.0.1
DOCKER_SWALO = $(DOCKER) agria-analytica/swalo:0.9.8-beta

DIRS = genome \
       results results/fastx results/fastqc \
       results/kmc results/kmc/tmp \
	   results/trimmomatic \
	   results/gatb-pipeline \
	   results/swalo \
	   results/quast \
	   results/ragtag results/ragtag/minia-k141_besst_ragtaq results/ragtag/minia-k141_ragtag results/minia-k141_swalo_ragtag
DATA = data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_1.fq.gz \
       data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_2.fq.gz

all: $(DIRS) qc trimmomatic genomescope gatb-pipeline swalo quast ragtag

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

trimmomatic: $(TRIM_FASTQ_PAIRED1)

$(TRIM_FASTQ_PAIRED1): $(DATA)
	$(DOCKER_TRIMMOMATIC) PE -threads 7 $^ $@ $(TRIM_FASTQ_UNPAIRED1) $(TRIM_FASTQ_PAIRED2) \
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
	$(DOCKER_KMC) kmc -k21 -t7 -m48 -ci3 -cs10000 $< $(basename $@) results/kmc/tmp/

# assemble genome.
gatb-pipeline: results/gatb-pipeline/kacang-lurik.assembly.fasta

# adjust working directory to avoid problem with minia unable to read/write files.
# to run multiple commands on docker we have to use /bin/bash -c or equivalent shell.
# source: https://stackoverflow.com/a/28490909
results/gatb-pipeline/kacang-lurik.assembly.fasta: $(TRIM_FASTQ_PAIRED1) $(TRIM_FASTQ_PAIRED2)
	$(DOCKER_GATB_PIPELINE) /bin/bash -c "cd $(dir $@) && ../../../gatb -1 ../../$(TRIM_FASTQ_PAIRED1) -2 ../../$(TRIM_FASTQ_PAIRED2) -o $(notdir $(basename $@)) --nb-cores 7" 

# scaffold using swalo.
ASSEMBLY_MINIA_K141 = results/gatb-pipeline/kacang-lurik.assembly_k141.contigs.fa
CONTIG_FILE = results/swalo/$(notdir $(ASSEMBLY_MINIA_K141))
BOWTIE2_IDX = $(CONTIG_FILE).1.bt2
MAP_FILE = $(basename $(CONTIG_FILE)).bowtie2.sam
MAX_INSERT_SIZE = 20000
MIN_CONTIG_LENGTH = 4400
swalo: results/swalo/scaffolds.fa

results/swalo/scaffolds.fa: $(CONTIG_FILE) results/swalo/unmappedOut.sam
	cd results/swalo && \
	$(DOCKER_SWALO) swalo $(notdir $<) $(MIN_CONTIG_LENGTH)

$(CONTIG_FILE): $(ASSEMBLY_MINIA_K141)
	ln -s $< $@

$(BOWTIE2_IDX): $(CONTIG_FILE) 
	$(DOCKER_SWALO) bowtie2-build --threads 4 $< $< 

$(MAP_FILE): $(TRIM_FASTQ_PAIRED1) $(TRIM_FASTQ_PAIRED2) $(BOWTIE2_IDX)
	$(DOCKER_SWALO) bowtie2 -k 5 -x $(CONTIG_FILE) -X $(MAX_INSERT_SIZE) -p 7 -1 $(TRIM_FASTQ_PAIRED1) -2 $(TRIM_FASTQ_PAIRED2) -S $@ 

results/swalo/myout.sam: $(MAP_FILE)
	cd results/swalo && \
	$(DOCKER_SWALO) bowtie2convert $(notdir $<) $(notdir $(CONTIG_FILE)) $(MAX_INSERT_SIZE)

results/swalo/unmappedOut.sam: $(CONTIG_FILE) results/swalo/myout.sam
	cd results/swalo && \
	$(DOCKER_SWALO) align $(notdir $<)

# assembly assessment using quast.
REF_AHYPOGAEA = genome/arahy.Tifrunner.gnm2.J5K5.genome_main.fna
REF_AHYPOGAEA_GFF = genome/arahy.Tifrunner.gnm2.ann1.4K0L.gene_models_main.gff3
ASSEMBLY = results/gatb-pipeline/kacang-lurik.assembly_k21.contigs.fa \
           results/gatb-pipeline/kacang-lurik.assembly_k41.contigs.fa \
		   results/gatb-pipeline/kacang-lurik.assembly_k61.contigs.fa \
		   results/gatb-pipeline/kacang-lurik.assembly_k81.contigs.fa \
		   results/gatb-pipeline/kacang-lurik.assembly_k101.contigs.fa \
		   results/gatb-pipeline/kacang-lurik.assembly_k121.contigs.fa \
		   results/gatb-pipeline/kacang-lurik.assembly_k141.contigs.fa \
		   results/gatb-pipeline/kacang-lurik.assembly.fasta
ASSEMBLY_NAME = k21,k41,k61,k81,k101,k121,k141,k141-scaffold
quast: results/gatb-pipeline/kacang-lurik.assembly.fasta
	$(DOCKER_QUAST) quast.py -o results/quast -t 8 --large -r $(REF_AHYPOGAEA) --features $(REF_AHYPOGAEA_GFF) -l $(ASSEMBLY_NAME) \
	--min-contig 500 $(ASSEMBLY)

# syntenic based scaffolding.
SCAFFOLD_BESST = results/gatb-pipeline/kacang-lurik.assembly.fasta
SCAFFOLD_OUTPUT = results/ragtag/minia-k141_besst_ragtag/ragtag.scaffold.fasta \
results/ragtag/minia-k141_ragtag/ragtag.scaffold.fasta
ragtag: $(SCAFFOLD_OUTPUT) 
	
results/ragtag/minia-k141_besst_ragtag/ragtag.scaffold.fasta: $(REF_AHYPOGAEA) $(SCAFFOLD_BESST)
	$(DOCKER_RAGTAG) ragtag.py scaffold -o $(dir $@) -t 7 -u $^

results/ragtag/minia-k141_ragtag/ragtag.scaffold.fasta: $(REF_AHYPOGAEA) $(ASSEMBLY_MINIA_K141)
	$(DOCKER_RAGTAG) ragtag.py scaffold -o $(dir $@) -t 7 -u $^



