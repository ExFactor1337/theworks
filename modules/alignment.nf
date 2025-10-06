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
        tuple val(ID), path("*bam"), emit: dedup_bam_ch
    script:
    opt_remove_duplicates = params.remove_duplicates ? '--REMOVE_DUPLICATES true' : ''
    """
    gatk MarkDuplicates --INPUT $BAM \
        --METRICS_FILE ${ID}.duplicate.metrics \
        --OUTPUT ${ID}.sorted.RG.MD.bam ${opt_remove_duplicates}
    """
}


process recalibrate_base_quality {
    tag "Recalibrating base quality scores for $ID"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(BAM), val(REF), path("*"), path("*") //"*" includes all reference files
        val(known_variants)

    output:
    // tuple val(meta), path(table)
    tuple val(ID), path("*recal_data.table"), emit: recal_table

    script:
    if (!known_variants) {
        throw new IllegalArgumentException("Error: The knownVariantsList for BaseRecalibrator cannot be empty.")
    }
    // Build the GATK command using the Groovy 'collect' and 'join' pattern
    def tabixCmds = known_variants.collect { path ->
        "tabix -f -p vcf \"${path}\""
    }.join('\n')
    def knownSitesArgs = known_variants.collect { path ->
        "--known-sites ${path}"
    }.join(' ')
    def recal_table = "${ID}.recal_data.table"
    """
    # Tabix the known variants
    ${tabixCmds}
    
    # Recalibrate base quality scores
    gatk BaseRecalibrator \
      --input ${BAM} \
      --reference ${REF}.fa \
      --output ${recal_table} ${knownSitesArgs}
    """
}

workflow{
    bwa_align(input_ch)
}
