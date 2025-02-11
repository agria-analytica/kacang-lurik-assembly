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
DOCKER_GENOMESCOPE = $(DOCKER) agria.analytica/genomescope2:be1953b

DIRS = genome \
	   data data/Fuhuasheng data/Shitouqi data/Tifrunner \
	   results results/fastx results/fastqc \
	   results/kmc results/kmc/tmp \
	   results/trimmomatic \
	   results/gatb-pipeline \
	   results/quast \
	   results/ragtag results/ragtag/minia-k141_besst_ragtaq results/ragtag/minia-k141_ragtag
	   results/genomescope
DATA = data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_1.fq.gz \
       data/X401SC21062291-Z01-F001/raw_data/KL/KL_DDSW210004672-1a_HCKFTDSX2_L1_2.fq.gz

all: $(DIRS) qc trimmomatic genomescope gatb-pipeline quast ragtag

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
	$(DOCKER_TRIMMOMATIC) PE -threads 10 $^ $@ $(TRIM_FASTQ_UNPAIRED1) $(TRIM_FASTQ_PAIRED2) $(TRIM_FASTQ_UNPAIRED2) \
	ILLUMINACLIP:$(ADAPTER):2:30:10 SLIDINGWINDOW:3:5

# genomescope
# k-mer table kacang lurik.
KMER_TABLE_21 = $(addsuffix .trim.k21.kmc_pre, $(addprefix results/kmc/, $(notdir $(basename $(DATA)))))
KMER_TABLE_31 = $(addsuffix .trim.k31.kmc_pre, $(addprefix results/kmc/, $(notdir $(basename $(DATA)))))
KMER_HIST_21  = results/kmc/kacang.lurik.k21.hist
KMER_HIST_31  = results/kmc/kacang.lurik.k31.hist
KMER_MERGE_21  = results/kmc/kacang.lurik.k21.merge.kmc_pre
KMER_MERGE_31  = results/kmc/kacang.lurik.k31.merge.kmc_pre

# k-mer table other A. hypogaea.
KMER_TABLE_TIFRUNNER_21 = results/kmc/Tifrunner.k21.kmc_pre
KMER_TABLE_TIFRUNNER_31 = results/kmc/Tifrunner.k31.kmc_pre
KMER_TABLE_FUHUASHENG_21 = results/kmc/Fuhuasheng.k21.kmc_pre
KMER_TABLE_FUHUASHENG_31 = results/kmc/Fuhuasheng.k31.kmc_pre
KMER_TABLE_SHITOUQI_21 = results/kmc/Shitouqi.k21.kmc_pre
KMER_TABLE_SHITOUQI_31 = results/kmc/Shitouqi.k31.kmc_pre
KMER_HIST_TIFRUNNER_21 = results/kmc/Tifrunner.k21.hist
KMER_HIST_TIFRUNNER_31 = results/kmc/Tifrunner.k31.hist
KMER_HIST_FUHUASHENG_21 = results/kmc/Fuhuasheng.k21.hist
KMER_HIST_FUHUASHENG_31 = results/kmc/Fuhuasheng.k31.hist
KMER_HIST_SHITOUQI_21 = results/kmc/Shitouqi.k21.hist
KMER_HIST_SHITOUQI_31 = results/kmc/Shitouqi.k31.hist

GENOMESCOPE_OUT= results/genomescope/kacang.lurik.k21.p2_linear_plot.png \
				 results/genomescope/kacang.lurik.k21.p3_linear_plot.png \
				 results/genomescope/kacang.lurik.k21.p4_linear_plot.png \
				 results/genomescope/kacang.lurik.k31.p2_linear_plot.png \
				 results/genomescope/kacang.lurik.k31.p3_linear_plot.png \
				 results/genomescope/kacang.lurik.k31.p4_linear_plot.png \
				 results/genomescope/Fuhuasheng.k21.p2_linear_plot.png \
				 results/genomescope/Fuhuasheng.k21.p4_linear_plot.png \
				 results/genomescope/Fuhuasheng.k31.p2_linear_plot.png \
				 results/genomescope/Fuhuasheng.k31.p4_linear_plot.png \
				 results/genomescope/Shitouqi.k21.p2_linear_plot.png \
				 results/genomescope/Shitouqi.k21.p4_linear_plot.png \
				 results/genomescope/Shitouqi.k31.p2_linear_plot.png \
				 results/genomescope/Shitouqi.k31.p4_linear_plot.png \
				 results/genomescope/Tifrunner.k21.p2_linear_plot.png \
				 results/genomescope/Tifrunner.k21.p4_linear_plot.png \
				 results/genomescope/Tifrunner.k31.p2_linear_plot.png \
				 results/genomescope/Tifrunner.k31.p4_linear_plot.png

genomescope: kmc $(GENOMESCOPE_OUT)
	
# run kmc.
kmc: $(KMER_HIST_21) $(KMER_HIST_31) $(KMER_HIST_FUHUASHENG_21) $(KMER_HIST_FUHUASHENG_31) $(KMER_HIST_SHITOQUI_21) $(KMER_HIST_SHITOQUI_31) $(KMER_HIST_TIFRUNNER_21) $(KMER_HIST_TIFRUNNER_31)

results/kmc/kacang.lurik.k%.hist: results/kmc/kacang.lurik.k%.merge.kmc_pre
	$(DOCKER_KMC) kmc_tools transform $(basename $<) histogram $@ -cx10000

results/kmc/%.k21.hist: results/kmc/%.k21.kmc_pre
	$(DOCKER_KMC) kmc_tools transform $(basename $<) histogram $@ -cx10000

results/kmc/%.k31.hist: results/kmc/%.k31.kmc_pre
	$(DOCKER_KMC) kmc_tools transform $(basename $<) histogram $@ -cx10000

# merge kmer from all fastq.
$(KMER_MERGE_21): $(KMER_TABLE_21)
	$(DOCKER_KMC) kmc_tools simple $(basename $^) union $(basename $@)

$(KMER_MERGE_31): $(KMER_TABLE_31)
	$(DOCKER_KMC) kmc_tools simple $(basename $^) union $(basename $@)

# generate kmer table from raw data.
$(KMER_TABLE_21): results/kmc/%.k21.kmc_pre: results/trimmomatic/%.gz
	$(DOCKER_KMC) kmc -k21 -t10 -m48 -ci3 -cs10000 $< $(basename $@) results/kmc/tmp/

$(KMER_TABLE_31): results/kmc/%.k31.kmc_pre: results/trimmomatic/%.gz
	$(DOCKER_KMC) kmc -k31 -t10 -m48 -ci3 -cs10000 $< $(basename $@) results/kmc/tmp/
	
results/kmc/Fuhuasheng.k21.kmc_pre: data/Fuhuasheng/Fuhuasheng.fa.gz
	$(DOCKER_KMC) kmc -k21 -t7 -m48 -ci3 -cs10000 -fa $< $(basename $@) results/kmc/tmp/

results/kmc/Fuhuasheng.k31.kmc_pre: data/Fuhuasheng/Fuhuasheng.fa.gz
	$(DOCKER_KMC) kmc -k31 -t7 -m48 -ci3 -cs10000 -fa $< $(basename $@) results/kmc/tmp/

results/kmc/Shitouqi.k21.kmc_pre: data/Shitouqi/Shitouqi.fa.gz
	$(DOCKER_KMC) kmc -k21 -t7 -m48 -ci3 -cs10000 -fa $< $(basename $@) results/kmc/tmp/

results/kmc/Shitouqi.k31.kmc_pre: data/Shitouqi/Shitouqi.fa.gz
	$(DOCKER_KMC) kmc -k31 -t7 -m48 -ci3 -cs10000 -fa $< $(basename $@) results/kmc/tmp/

results/kmc/Tifrunner.k21.kmc_pre: data/Tifrunner/Tifrunner.fa.gz
	$(DOCKER_KMC) kmc -k21 -t7 -m48 -ci3 -cs10000 -sm -fa $< $(basename $@) results/kmc/tmp/

results/kmc/Tifrunner.k31.kmc_pre: data/Tifrunner/Tifrunner.fa.gz
	$(DOCKER_KMC) kmc -k31 -t7 -m48 -ci3 -cs10000 -sm -fa $< $(basename $@) results/kmc/tmp/

# download data from ENA using wget, strip their quality values, save as compressed fasta.
data/Fuhuasheng/Fuhuasheng.fa.gz: data/Fuhuasheng/fastq.list
	wget -O - -i $< | gunzip -c | perl -lane 'if (($$.%4)==1) {print ">".substr($$_,1)};if(($$.%4)==2){print $$_}' | pigz -p8 -c > $@

data/Shitouqi/Shitouqi.fa.gz: data/Shitouqi/fastq.list
	wget -O - -i $< | gunzip -c | perl -lane 'if (($$.%4)==1) {print ">".substr($$_,1)};if(($$.%4)==2){print $$_}' | pigz -p8 -c > $@

data/Tifrunner/Tifrunner.fa.gz: data/Tifrunner/fastq.list
	wget -O - -i $< | gunzip -c | perl -lane 'if (($$.%4)==1) {print ">".substr($$_,1)};if(($$.%4)==2){print $$_}' | pigz -p8 -c > $@

# run genomescope.
results/genomescope/kacang.lurik.k%.p2_linear_plot.png: results/kmc/kacang.lurik.k%.
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 2 --kmercov 5

results/genomescope/kacang.lurik.k%.p3_linear_plot.png: results/kmc/kacang.lurik.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 3 --kmercov 5

results/genomescope/kacang.lurik.k%.p4_linear_plot.png: results/kmc/kacang.lurik.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 4 --kmercov 5

results/genomescope/Fuhuasheng.k%.p2_linear_plot.png: results/kmc/Fuhuasheng.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 2 --kmercov 5

results/genomescope/Fuhuasheng.k%.p4_linear_plot.png: results/kmc/Fuhuasheng.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 4 --kmercov 5

results/genomescope/Shitouqi.k%.p2_linear_plot.png: results/kmc/Shitouqi.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 2 --kmercov 25

results/genomescope/Shitouqi.k%.p4_linear_plot.png: results/kmc/Shitouqi.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 4 --kmercov 25

results/genomescope/Tifrunner.k%.p2_linear_plot.png: results/kmc/Tifrunner.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 2 --kmercov 50

results/genomescope/Tifrunner.k%.p4_linear_plot.png: results/kmc/Tifrunner.k%.hist
	$(DOCKER_GENOMESCOPE) genomescope2 -i $< -o $(dir $@) -k $* -n $(firstword $(subst _, ,$(notdir $@))) -p 4 --kmercov 50

# assemble genome.
gatb-pipeline: results/gatb-pipeline/kacang-lurik.assembly.fasta

# adjust working directory to avoid problem with minia unable to read/write files.
# to run multiple commands on docker we have to use /bin/bash -c or equivalent shell.
# source: https://stackoverflow.com/a/28490909
results/gatb-pipeline/kacang-lurik.assembly.fasta: $(TRIM_FASTQ_PAIRED1) $(TRIM_FASTQ_PAIRED2)
	$(DOCKER_GATB_PIPELINE) /bin/bash -c "cd $(dir $@) && ../../../gatb -1 ../../$(TRIM_FASTQ_PAIRED1) -2 ../../$(TRIM_FASTQ_PAIRED2) -o $(notdir $(basename $@)) --nb-cores 7" 

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
ASSEMBLY_MINIA_K141 = results/gatb-pipeline/kacang-lurik.assembly_k141.contigs.fa
SCAFFOLD_OUTPUT = results/ragtag/minia-k141_besst_ragtag/ragtag.scaffold.fasta \
results/ragtag/minia-k141_ragtag/ragtag.scaffold.fasta
ragtag: $(SCAFFOLD_OUTPUT) 
	
results/ragtag/minia-k141_besst_ragtag/ragtag.scaffold.fasta: $(REF_AHYPOGAEA) $(SCAFFOLD_BESST)
	$(DOCKER_RAGTAG) ragtag.py scaffold -o $(dir $@) -t 7 -u $^

results/ragtag/minia-k141_ragtag/ragtag.scaffold.fasta: $(REF_AHYPOGAEA) $(ASSEMBLY_MINIA_K141)
	$(DOCKER_RAGTAG) ragtag.py scaffold -o $(dir $@) -t 7 -u $^


