# NCCR Genomics Pipeline

**Under Development**


Pipeline for analysis of isolate genomes. 

### Steps Currently Included:

- QC
- Isolate Genome Assembly 
- Gene calling and annotation
- Mapping back to assembly or to reference genome
- Variant calling and annotation
- ANI calculation
- Phylogenetics with PhyloPhlAn
- Pangenome analysis with panX


### Installation

The software depends on `conda`, `Python >= 3.8` and `snakemake == 5.22`. 

```
git clone https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline
cd NCCR_genomicsPipeline/code/package
conda env create -f nccrPipe_environment.yaml
# mamba env create -f nccrPipe_environment.yaml
conda activate nccrPipe
pip install -e . 

```

### Running the pipeline

```
nccrPipe COMMAND [OPTIONS]

Commands:
  isolate
  rnaseq
  unlock

Options:
  -c, --config TEXT Path to configuration file 
  -m, --method TEXT Method/Analysis to run
  --dry             Show commands without running them
  --local           Run locally/without submitting to the cluster
  --help            Show the help message

```

### Running RNAseq pipeline

- To run STAR/featureCounts pipeline run:

```
nccrPipe rnaseq -c <configfile> -m star

```

To see what jobs are going to be submitted to the cluster add `--dry` flag. To run pipeline without submitting jobs to the cluster add `--local`` flag.


- To run kallisto pipeline run:

```

nccrPipe rnaseq -c <configfile> -m kallisto

```

### Running isolate pipeline

UNDER CONSTRUCTION


### RNAseq Config File:

- YAML file with following mandatory fields:

```
ProjectName: Test
dataDir: path to data directory (see data structure below)
outDir: path to output directory
sampleFile: file with sample names (see example below)
fq_fwd: _1.fq.gz # forward reads fastq suffix
fq_rvr: _2.fq.gz # reverse reads fastq suffix

#Preprocessing
qc: yes
mink: 11
trimq: 14 # 37
mapq: 20
minlen: 45

merged: false # By default to do not use merge for isolate genome assembly
fastqc: no # Options: no, before, after, both


# STAR
refGenome: reference_genome.fna # Uncompressed (Did not test with compressed file)
refAnn: reference_annotation.gtf # Important to have .gtf not a .gff
genomeDir: directory to output genome index
overhang: 149 # ReadLength - 1
maxIntron: 50000 # Max size of intron (Depends on organism) 


# featureCounts
strand: 0 # Strandiness of the RNASeq, can be 0,1,2
attribute: gene_id # Should work, if using .gtf file 


# kallisto
transcriptome: transcriptome.fa #Uncompressed
kallistoIdx: transcriptome_index_file.idx


# Standard parameters. Dont change these!!!
adapters: '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/resources/adapters.fa'
phix: '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/resources/phix174_ill.ref.fa.gz'
bbmap_human_ref: '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/resources/human_bbmap_ref/'
bbmap_mouse_ref: '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/resources/mouse_bbmap_ref/'

```

### Example Data Structure 

```
|--data/
     |--samples.txt   
     |--raw/
     |   |--Sample1
     |   |      |--Sample1_abcd_R1.fq.gz
     |   |      |--Sample1_abcd_R2.fq.gz
     |   |
     |   |--Sample2
     |          |--Sample1_efgh_R1.fq.gz
     |          |--Sample2_efgh_R2.fq.gz 
     |   
     |--processed/     

```
### Example samples.txt:
```
Sample1
Sample2
```

### Example config

```
ProjectName: Example Config
dataDir: data/raw
outDir: data/processed
sampleFile: data/samples.txt
fq_fwd: _R1.fq.gz 
fq_rvr: _R2.fq.gz 

...

```
- Test RNASeq config is `code/package/nccrPipe/configs/rnaseq_config.yaml`
- Test samples file is `code/package/nccrPipe/configs/rnaseq_config.yaml`
- RNAseq snakemake rules are in `code/package/nccrPipe/rnaseq_rules`

----------------------------------------------------------
# UNDER CONSTRUCTION 

### Analysis Options:

- preprocess
- merge_fastq
- assemble
- quastCheck
- plasmid
- annotate
- nucmer
- runANI
- phylogeny
- findAb
- pileup
- align
- align_with_ref
- call_vars
- call_vars2 -> use this one
- call_vars3
- markdup
- type
- serotype
- anVar
- anVar2
- pangenome

All of these need to be tested and documented

## Example

Default data structure: 


By default, data will be in `data/raw` and the output directory `data/processed`

```
cd project
nccrPipe -a create_config -c code/configs/project_config.yaml

```


## Variant Calling

### General Steps:
1. Align reads to reference genome with BWA. 
2. Remove duplicates with GATK MarkDuplicates
3. Run `bcftools mpileup` + `bcftools call`
4. `bcftools filter` 

```
-g5    filter SNPs within 5 base pairs of an indel or other other variant type
-G10   filter clusters of indels separated by 10 or fewer base pairs allowing only one to pass
-e  exclude
QUAL<10 calls with quality score < 10
DP4[2]<10 || DP4[3]<10  calls with < 10 reads (forward and reverse) covering the variant
(DP4[2] + DP4[3])/sum(DP4) < 0.9 calls with allele frequence < 90 %
MQ<50 calls with average mapping quality < 50 

```  
- Add filter based on coverage? Regions with really high coverage generally contain lots of artifacts. 

### Filters:


### Annotation:

- Still a little complicated. Adding new genome to the snpEff database not integrated into the main pipeline.
- Created conda environment within the pipeline: `.snakemake/conda/`
- If want to add new genome, run
```
nccrPipe -c <config_file> -a addGenomeSnpEff
```

- Config file has to include path to the gbk (and ideally fasta) files, as well as the name of the genome. The chromosome names between the files have to match. 
- If that doesn't fail, can run

```
nccrPipe -c <config_file> -a anVar
```

a
