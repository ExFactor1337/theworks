nextflow.enable.dsl=2

process run_fastp {
    tag "Running fastp on $ID"
    publishDir params.out_dir, mode: 'symlink'
    input: 
        tuple val(ID), path(FQ1), path(FQ2)
    output:
        tuple val(ID), path("*trimmed_1.fastq"), path("*trimmed_2.fastq"), emit: trim_fastq_ch
        path("*html")
        path("*json")
    script:
    """
    fastp --in1 ${FQ1} \
          --in2 ${FQ2} \
          --out1 ${ID}.trimmed_1.fastq \
          --out2 ${ID}.trimmed_2.fastq \
          --thread ${params.fastp_cpus} \
          --html ${ID}_fastp_report.html \
          --json ${ID}_fastp_report.json
    """
}

