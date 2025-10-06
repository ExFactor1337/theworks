nextflow.enable.dsl=2
import Utils
import ParamsChecker
import Auditor

include { run_fastp } from './modules/read_qc.nf'
include { create_sequence_dictionary; 
          bwa_align; 
          add_read_groups; 
          mark_duplicates;
          recalibrate_base_quality } from './modules/alignment.nf'

fastq_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2}.{fastq,fq}{.gz,}", flat: true)

//Directory of bwa-indexed reference files 
Utils.validateReferenceFiles(params.reference_dir)
ref_ch = Channel.fromPath("${params.reference_dir}/*")
            .collect()
            .map { tuple Utils.getCommonPrefix(it), it }

// Collect all files in the known_variants directory
known_variants_ch = Channel.fromPath("${params.known_variants_dir}/*{vcf,vcf.gz}")
                            .collect()


// Pass in the params 
def auditor = new Auditor()

workflow {
    // Handle generation of gatk sequence dictionary if not provided
    // and update the ref_ch to include the dict file
    if(!params.gatk_sequence_dict) {
        create_sequence_dictionary(ref_ch)
        ref_ch = ref_ch.combine(create_sequence_dictionary.out.dict_ch).map { id, file_list, dict_path ->
                def updated_file_list = file_list + [dict_path]
                return tuple(id, updated_file_list)
        }
    }
    
    // Quality check fastqs
    run_fastp(fastq_ch)

    // Align fastqs to reference
    bwa_align(
        run_fastp.out.trim_fastq_ch
            .combine(ref_ch)
    )
    
    // Add read groups
    add_read_groups(bwa_align.out.sorted_bam_ch)

    // Mark duplicates
    mark_duplicates(add_read_groups.out.rg_bam_ch)

    // Recalibrate base quality scores
    recalibrate_base_quality(mark_duplicates.out.dedup_bam_ch
                                .combine(ref_ch)
                                .combine(known_variants_ch),
                             known_variants_ch)
}

workflow.onComplete {
    auditor.generateAuditLog(
        workflow, 
        params.findAll { k, v -> !(v instanceof Map) }, 
        params.ascii_art
    )
}
