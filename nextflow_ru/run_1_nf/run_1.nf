#!/usr/bin/env nextflow

// 20230113 I am rerunning this pipeline with an alternative transcriptome at Ru's request.
// See the comment in run_1_nf.r for more details.

/*
A Nextflow pipeline for processing the run1 data of Ru.
It is single-end bulk RNA-seq
*/

samples_ch = Channel.fromFilePairs("/home/humebc/projects/ru/raw_seq_data/run_1/*.fastq.gz", size: 1).map{[it[0], it[1][0]]}
reference = file("/home/humebc/projects/ru/reference/alt_transcriptome/Phaeodactylum_tricornutum.ASM15095v2.cds.all.kallisto.index")

process fastp{
    tag "${sample}"
    conda "fastp"
    publishDir "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results_alt_ref/fastp", pattern: "*.html"

    input:
    tuple val(sample), path(read_1) from samples_ch

    output:
    file "${sample}.fastp.html" into fastp_out_ch
    tuple val(sample), file("${sample}.clean.fq.gz") into kallisto_in_ch

    script:
    """
    fastp -q 20 -i $read_1 -o ${sample}.clean.fq.gz
    mv fastp.html ${sample}.fastp.html
    """
}

process kallisto{
    tag "${sample}"
    container "jennylsmith/kallistov45.0:latest"
    publishDir "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results_alt_ref/kallisto/${sample}"
    cpus 10

    input:
    tuple val(sample), path(read_1_clean) from kallisto_in_ch
    file reference

    output:
    file "*" into kallisto_out_ch

    script:
    """
    kallisto quant -i $reference -o . -t ${task.cpus} --single -l 200 -s 30 $read_1_clean
    """
}
