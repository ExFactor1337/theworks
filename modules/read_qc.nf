nextflow.enable.dsl=2

params.fastq_dir = false
params.out_dir = 'results'

fastq_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2,1.clipped,2.clipped}.{fastq,fq}{.gz,}", flat: true)


process run_fastp {
    publishDir params.out_dir, mode: 'symlink'
    input: 
        tuple val(ID), path(FQ1), path(FQ2)
    output:
        path("*")
    script:
    """
    fastp --in1 ${FQ1} \
          --in2 ${FQ2} \
          --out1 ${ID}.trimmed_1.fastq \
          --out2 ${ID}.trimmed_2.fastq \
          --thread 6 \
          --html ${ID}_fastp_report.html \
          --json ${ID}_fastp_report.json
    """
}

workflow {
    run_fastp(fastq_ch)
}