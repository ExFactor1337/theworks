nextflow.enable.dsl=2

import Utils
include { run_fastp } from './modules/read_qc.nf'
include { bwa_align } from './modules/alignment.nf'

//fastq_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2,1.clipped,2.clipped}.{fastq,fq}{.gz,}", flat: true)

println Utils.summarizeRun(workflow, params, log)

workflow {
    //run_fastp(fastq_ch)
    //bwa_align()
}


