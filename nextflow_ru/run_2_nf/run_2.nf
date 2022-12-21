#!/usr/bin/env nextflow

/*
A Nextflow pipeline for processing the run1 data of Ru.
It is single-end bulk RNA-seq
*/

samples_ch = Channel.fromFilePairs("/home/humebc/projects/ru/raw_seq_data/run_2/X204SC22093327-Z01-F001_0*/01.RawData/*/*{_1,_2}.fq.gz", size: 2)

reference = file("/home/humebc/projects/ru/reference/GCF_000150955.2_ASM15095v2_rna.kallisto.index")

process fastp{
    tag "${sample}"
    conda "fastp"
    publishDir "/home/humebc/projects/ru/nextflow_ru/run_2_nf/results/fastp", pattern: "*.html"

    input:
    tuple val(sample), path(reads) from samples_ch

    output:
    file "${sample}.fastp.html" into fastp_out_ch
    tuple val(sample), file("${sample}.clean.1.fq.gz"), file("${sample}.clean.2.fq.gz")  into kallisto_in_ch

    script:
    """
    fastp -q 20 -i ${reads[0]} -I ${reads[1]} -o ${sample}.clean.1.fq.gz -O ${sample}.clean.2.fq.gz
    mv fastp.html ${sample}.fastp.html
    """
}

process kallisto{
    tag "${sample}"
    container "jennylsmith/kallistov45.0:latest"
    publishDir "/home/humebc/projects/ru/nextflow_ru/run_2_nf/results/kallisto/${sample}"
    cpus 10

    input:
    tuple val(sample), path(read_1_clean), path(read_2_clean) from kallisto_in_ch
    file reference

    output:
    file "*" into kallisto_out_ch

    script:
    """
    kallisto quant -i $reference -o . -t ${task.cpus} $read_1_clean $read_2_clean
    """
}