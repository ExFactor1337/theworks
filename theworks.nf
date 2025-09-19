nextflow.enable.dsl=2
import Utils
import ParamsCheck

include { run_fastp } from './modules/read_qc.nf'
include { bwa_align; add_read_groups; mark_duplicates } from './modules/alignment.nf'


fastq_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2,1.clipped,2.clipped}.{fastq,fq}{.gz,}", flat: true)

//Directory of bwa-indexed reference files 
ref_ch = Channel.fromPath("${params.reference_dir}/*")
            .collect()
            .map { tuple Utils.getCommonPrefix(it), it }

//println Utils.summarizeRun(workflow, params, log)

//println params

// Instantiate ParamsCheck and create vparams object to hold its output
def paramsChecker = new ParamsCheck(params, params.definitions)
def vparams = paramsChecker.validateAndBuild()

println vparams