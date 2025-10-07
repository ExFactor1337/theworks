nextflow.enable.dsl=2
process call_variants_haplotypecaller {
    tag "Calling variants for $ID"
    publishDir params.out_dir, mode: 'symlink'
    input:
    	tuple val(ID), path(BAM), path(BAI), val(REF), path("*")

    output:
    	tuple val(ID), path("${ID}.gvcf.gz"), path("${ID}.gvcf.gz.tbi"), emit: raw_vcf_ch

    script:
    def output_gvcf = "${ID}.gvcf.gz"
    """
    gatk HaplotypeCaller \
      --reference ${REF}.fa \
      --input ${BAM} \
      --output ${output_gvcf} \
      --emit-ref-confidence GVCF \
      --sample-name ${ID}
    """
}

process generate_snp_recalibration_model {
    tag "Generating SNP recalibration model for $ID"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(VCF), val(REF), path("*"), path("*") //"*" includes all reference files

    output:
        tuple val(ID), path("*snp.recal"), path("*snp.tranches"), path("*snp.plots.R"), emit: snp_recal_ch

    script:
    def recal_file = "${ID}.snp.recal"
    def tranches_file = "${ID}.snp.tranches"
    def plots_r_file = "${ID}.snp.plots.R"
    """
    gatk VariantRecalibrator \
      --input ${VCF} \
      --reference ${REF}.fa \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
      --resource:omni,known=false,training=true,truth=true,prior=12.0 ${params.omni} \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.thousand_genomes} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
      -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an MQ \
      -mode SNP \
      --tranches-file ${tranches_file} \
      --rscript-file ${plots_r_file} \
      --output ${recal_file}
    """
}

process generate_indel_recalibration_model {
    tag "Generating INDEL recalibration model for $ID"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(VCF), val(REF), path("*"), path("*") //"*" includes all reference files

    output:
        tuple val(ID), path("*indel.recal"), path("*indel.tranches"), path("*indel.plots.R"), emit: indel_recal_ch

    script:
    def recal_file = "${ID}.indel.recal"
    def tranches_file = "${ID}.indel.tranches"
    def plots_r_file = "${ID}.indel.plots.R"
    """
    gatk VariantRecalibrator \
      --input ${VCF} \
      --reference ${REF}.fa \
      // INTEGRATE KNOWN INDEL RESOURCES AS NEEDED
      -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an MQ \
      -mode INDEL \
      --tranches-file ${tranches_file} \
      --rscript-file ${plots_r_file} \
      --output ${recal_file}
    """
} 

process apply_snp_recalibration {
    tag "Applying SNP recalibration for $ID"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(VCF), path(VCF_IDX), path(SNP_RECAL), path(SNP_TRANCHES), val(REF), path("*"), path("*") //"*" includes all reference files

    output:
        tuple val(ID), path("${ID}.snp.recal.vcf.gz"), path("${ID}.snp.recal.vcf.gz.tbi"), emit: snp_filtered_vcf_ch

    script:
    def output_vcf = "${ID}.snp.recal.vcf.gz"
    """
    gatk ApplyVQSR \
      --input ${VCF} \
      --reference ${REF}.fa \
      --truth-sensitivity-filter-level 99.5 \
      --tranches-file ${SNP_TRANCHES} \
      --recal-file ${SNP_RECAL} \
      -mode SNP \
      --output ${output_vcf}
    """
} 

process apply_indel_recalibration {
    tag "Applying INDEL recalibration for $ID"
    publishDir params.out_dir, mode: 'symlink'
    input:
        tuple val(ID), path(VCF), path(VCF_IDX), path(INDEL_RECAL), path(INDEL_TRANCHES), val(REF), path("*"), path("*") //"*" includes all reference files

    output:
        tuple val(ID), path("${ID}.indel.recal.vcf.gz"), path("${ID}.indel.recal.vcf.gz.tbi"), emit: final_vcf_ch

    script:
    def output_vcf = "${ID}.indel.recal.vcf.gz"
    """
    gatk ApplyVQSR \
      --input ${VCF} \
      --reference ${REF}.fa \
      --truth-sensitivity-filter-level 99.0 \
      --tranches-file ${INDEL_TRANCHES} \
      --recal-file ${INDEL_RECAL} \
      -mode INDEL \
      --output ${output_vcf}
    """
} 