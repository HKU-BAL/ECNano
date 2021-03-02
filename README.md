# ECNano
A target captured medical exome ONT sequencing with amplicon variant calling workflow

## Introduction

Having a high sequencing error of ONT and limited throughput from a single MinION flowcell ultimately limits its applicability for accurate variant detection. Medical exome sequencing (MES) targets clinically significant exon regions, allowing rapid and comprehensive screening of pathogenic variants. By applying MES with MinION sequencing, the technology can achieve a more uniform capture of the target regions, shorter turnaround time, and lower sequencing cost per sample. 

ECNano is an out-of-the-box workflow comprising (1) a wet-lab protocol for ONT target enrichment sequencing and (2) a downstream variant detection workflow, Clair-ensemble. The ECNano wet-lab protocol was optimized to perform long-read target enrichment and ONT library preparation to stably generate high-quality MES data with adequate coverage. The subsequent variant-calling workflow, Clair-ensemble, adopted a fast RNN-based variant caller, Clair, and was optimized for target enrichment data. ECNano was evaluated on both reference DNA samples and patient samples.

## Build an anaconda virtual environment
Please install anaconda using the installation guide at https://docs.anaconda.com/anaconda/install/
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the clair-ensemble environment and install dependencies
conda create -n clair-ensemble -c bioconda -y clair
conda activate clair-ensemble
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree
conda install porechop minimap2 samtools snakemake bioawk

git clone --depth 1 https://github.com/HKU-BAL/ECNano

# unset the perl env variables pointing to different perl installations
unset PERL5LIB PERL_LOCAL_LIB_ROOT 

# install guppy if input fast5
cd ECNano
curl -O https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.4.4_linux64.tar.gz && tar -xf ont-guppy_3.4.4_linux64.tar.gz
```
## Workflow usage
```
snakemake --cores INT [--config discard_middle=True|False 
                                minionqc=True|False 
                                gpu_device=INT 
                                in_dir=DIR 
                                out_dir=DIR
                                BED_FILE_PATH=PATH
                                REFERENCE_FASTA_FILE_PATH=PATH]
```
For all available configs, please refer to `config.yaml`

## Demo run
```
mkdir input reference
curl -o reference/GRCh38_no_alt_analysis_set.no_chr.fasta http://www.bio8.cs.hku.hk/dataset/ECNano/GRCh38_no_alt_analysis_set.no_chr.fasta \
        input/HG001_ECNano_dataset1.fastq http://www.bio8.cs.hku.hk/dataset/ECNano/HG001_ECNano_dataset1.fastq
samtools faidx reference/GRCh38_no_alt_analysis_set.no_chr.fasta
snakemake --cores 24
```
