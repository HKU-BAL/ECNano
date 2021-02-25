import os
from glob import glob
configfile: "config.yaml"

rule all:
    input:
        expand("{out_dir}/Clair_ensemble_out",out_dir=config["out_dir"])
fastq_list=glob("input/*.f*q")
if fastq_list:
    rule fastq_in:
        input:
            "input"
        output:
            fq=expand("{out_dir}/all.pass.fq",out_dir=config["out_dir"])
        shell:
            """bioawk -c fastx '{{if (meanqual($qual)>{config[qual_threshold]}) print "@"$name" "$comment"\\n"$seq"\\n+\\n"$qual}}' {input}/*.f*q > {output.fq} """

fast5_list=glob("input/*.f*5")
if fast5_list:
    rule fast5_in:
        input: 
            "input"
        output:
            expand("{out_dir}/all.pass.fq",out_dir=config["out_dir"])
        shell:
            "bin/ont-guppy/bin/guppy_basecaller -r -i {input} -s {config[out_dir]}/Guppy_out -c dna_r9.4.1_450bps_hac.cfg --device cuda:{config[gpu_device]} -q 0 --hp_correct 1 --qscore_filtering --min_qscore 3 && "
            "cat {config[out_dir]}/Guppy_out/pass/*.fastq > {config[out_dir]}/all.pass.fq"

rule analysis:
    input:
        expand("{out_dir}/all.pass.fq",out_dir=config["out_dir"])
    output:
        touch(".analysis.done")
    shell:
        """echo "[INFO] *average read quality: " $(bioawk -c fastx '{{ print meanqual($qual) }}' {input} | awk '{{sum += $1}} END {{print sum/NR}}');"""
        """echo "[INFO] *total throughput: " $(cat {input} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c);"""
        """echo "[INFO] *average read length: " $(cat {input} | awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}')"""

if config["minionqc"]==True:
    rule minionqc:
        input:
            expand("{out_dir}/Guppy_out/sequencing_summary.txt",out_dir=config["out_dir"])
        output:
            expand("{out_dir}/MinIONQC_out",out_dir=config["out_dir"])
        run:
            directory(shell("MinIONQC.R -i {input} -o {output} -p {config[threads]}"))

rule porechop:
    input:
        checkpoint=".analysis.done",
        fq=expand("{out_dir}/all.pass.fq",out_dir=config["out_dir"])
    output:
        expand("{out_dir}/Porechop_out/guppy_pass.porechop.fastq",out_dir=config["out_dir"])
    run:
        if config["discard_middle"]==True:
            discard_middle_flag="--discard_middle"
        else:
            discard_middle_flag=""   
        shell("mkdir -p {config[out_dir]}/Porechop_out && " \
        "porechop -i {input.fq} -o {output} -t {config[threads]} {discard_middle_flag}")

rule align_index:
    input:
        expand("{out_dir}/Porechop_out/guppy_pass.porechop.fastq",out_dir=config["out_dir"])
    output:
        bam=expand("{out_dir}/Minimap2_out/guppy_pass.porechop.minimap2_hg38.sorted.bam",out_dir=config["out_dir"]),
        bai=expand("{out_dir}/Minimap2_out/guppy_pass.porechop.minimap2_hg38.sorted.bam.bai",out_dir=config["out_dir"])
    shell:
        "mkdir -p {config[out_dir]}/Minimap2_out && "
        "minimap2 -ax map-ont reference/GRCh38_no_alt_analysis_set.no_chr.fasta {input} -t {config[threads]} | samtools sort -o {output.bam}; "
        "samtools index {output.bam}"

rule clair_ensemble:
    input:
        bam=expand("{out_dir}/Minimap2_out/guppy_pass.porechop.minimap2_hg38.sorted.bam",out_dir=config["out_dir"]),
        bai=expand("{out_dir}/Minimap2_out/guppy_pass.porechop.minimap2_hg38.sorted.bam.bai",out_dir=config["out_dir"])
    output:
        directory(expand("{out_dir}/Clair_ensemble_out",out_dir=config["out_dir"]))
    run:
        config["CLAIR_MODELS"]= [''.join([os.getcwd(),'/',path]) for path in config["CLAIR_MODELS"]]
        clair_models_delim=','.join(config["CLAIR_MODELS"])
        shell("./bin/run_Clair_ensemble.sh -b {input.bam} -d {config[BED_FILE_PATH]} -r {config[REFERENCE_FASTA_FILE_PATH]} -c {config[CLAIR]} -m {clair_models_delim} -e {config[ENSEMBLE_CPP_EXECUTABLE]} -t {config[threads]} -o {output}")




