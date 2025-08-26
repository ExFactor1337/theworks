nextflow.enable.dsl=2

process bwa_align {
    tag "Aligning $ID against $REF"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(FQ1), path(FQ2), val(REF), path("*")
    output:
        tuple val(ID), path("*sorted.bam"), emit: orig_bam_ch
    script:
    """
        bwa mem -t 4 $REF $FQ1 $FQ2 | \
            samtools sort -o "${ID}_${REF}.sorted.bam"
    """
}

process add_read_groups {
    tag "Adding read group tags to $ID"
    publishDir params.out_dir, mode: 'symlink'
    input: 
        tuple val(ID), path(BAM)
    output: 
        tuple val(ID), path("*sortedRG.bam"), emit: rg_bam_ch
    script:
    """
    # LB, PL, etc are required inputs that need dummying in. 
    gatk AddOrReplaceReadGroups -I $BAM -O ${ID}.sortedRG.bam \
        -LB . \
        -PL . \
        -PU . \
        -SM ${ID}
    """
}

process mark_duplicates {
    input:
        tuple val(ID), path(BAM) 
    output:
        tuple val(ID), path
    script:
    """
    gatk MarkDuplicates -I 
    """
}

workflow{
    bwa_align(input_ch)
}
