nextflow.enable.dsl=2
import Utils
import ParamsChecker
import Auditor

include { run_fastp } from './modules/read_qc.nf'
include { create_sequence_dictionary; bwa_align; add_read_groups; mark_duplicates } from './modules/alignment.nf'

fastq_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2}.{fastq,fq}{.gz,}", flat: true)

//Directory of bwa-indexed reference files 
ref_ch = Channel.fromPath("${params.reference_dir}/*")
            .collect()
            .map { tuple Utils.getCommonPrefix(it), it }.view()

// Pass in the params 
def auditor = new Auditor(workflow, params.findAll { k, v -> !(v instanceof Map) })
/*
workflow {
    if(!params.gatk_sequence_dict) {
        create_sequence_dictionary(ref_ch)
        ref_ch = ref_ch.combine(create_sequence_dictionary.out.dict_ch).map { id, file_list, dict_path ->
                def updated_file_list = file_list + [dict_path]
                return tuple(id, updated_file_list)
        }
    }
    
    /*
    // Quality check fastqs
    run_fastp(fastq_ch)

    // Align fastqs to reference
    bwa_align(
        run_fastp.out.trim_fastq_ch
            .combine(ref_ch)
    )
    
}
*/