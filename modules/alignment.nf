nextflow.enable.dsl=2

process bwa_align {
    tag "Aligning $ID against $REF"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(FQ1), path(FQ2), val(REF), path("*")
    output:
        path("*")
    script:
    """
        bwa mem -t 4 $REF $FQ1 $FQ2 | \
            samtools sort -o "${ID}_${REF}.sorted.bam"
        samtools index "${ID}_${REF}.sorted.bam"
    """
}

workflow{
    bwa_align(input_ch)
}
