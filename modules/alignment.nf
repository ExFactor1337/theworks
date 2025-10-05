nextflow.enable.dsl=2

process create_sequence_dictionary {
    tag "Creating sequence dictionary for $REF"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(REF), path("*")
    output:
        path("*dict"), emit: dict_ch
    script:
    """
    gatk CreateSequenceDictionary -R ${REF}.fa -O ${REF}.dict
    """
}   

process bwa_align {
    tag "Aligning $ID against $REF"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(FQ1), path(FQ2), val(REF), path("*")
    output:
        tuple val(ID), path("*sorted.bam"), emit: sorted_bam_ch
    script:
    """
        bwa mem -t ${params.bwa_cpus} $REF $FQ1 $FQ2 | \
            samtools sort -o "${ID}_${REF}.sorted.bam"
    """
}

process add_read_groups {
    tag "Adding read group tags to $ID"
    publishDir params.out_dir, mode: 'symlink'
    input: 
        tuple val(ID), path(BAM)
    output: 
        tuple val(ID), path("*sorted.RG.bam"), emit: rg_bam_ch
    script:
    """
    # LB, PL, etc are required inputs that need dummying in. 
    gatk AddOrReplaceReadGroups -I $BAM -O ${ID}.sorted.RG.bam \
        -LB . \
        -PL . \
        -PU . \
        -SM ${ID}
    """
}

process mark_duplicates {
    tag "Marking duplicates in $ID"
    publishDir params.out_dir, mode: 'symlink' 
    input:
        tuple val(ID), path(BAM) 
    output:
        tuple val(ID), path("*bam"), path("*metrics"), emit: dedup_bam_ch
    script:
    opt_remove_duplicates = params.remove_duplicates ? '--REMOVE_DUPLICATES true' : ''
    """
    gatk MarkDuplicates --INPUT $BAM \
        --METRICS_FILE ${ID}.duplicate.metrics \
        --OUTPUT ${ID}.sorted.RG.MD.bam ${opt_remove_duplicates}
    """
}

// Common parameters for BQSR
params.fasta = 'path/to/reference.fasta'
params.known_snps = 'path/to/dbsnp.vcf.gz'
params.known_indels = 'path/to/mills_indels.vcf.gz'

process recalibrate_base_quality {
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(BAM), val(REF), path("*") //"*" includes all reference files
        path(known_snps)
        path(known_indels)

    output:
    // tuple val(meta), path(table)
    tuple val(ID, path("*recal_data.table")), emit: recal_table

    script:
    def recal_table = "${meta.id}.recal_data.table"
    """
    gatk BaseRecalibrator \\
      --input ${input_bam} \\
      --reference ${ref} \\
      --output ${recal_table} \\
      --known-sites ${known_snps} \\
      --known-sites ${known_indels}
    """
}

workflow{
    bwa_align(input_ch)
}
