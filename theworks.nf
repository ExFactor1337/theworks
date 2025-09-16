nextflow.enable.dsl=2
import Utils

include { run_fastp } from './modules/read_qc.nf'
include { bwa_align; add_read_groups; mark_duplicates } from './modules/alignment.nf'


fastq_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2,1.clipped,2.clipped}.{fastq,fq}{.gz,}", flat: true)

//Directory of bwa-indexed reference files 
ref_ch = Channel.fromPath("${params.reference_dir}/*")
            .collect()
            .map { tuple Utils.getCommonPrefix(it), it }

//println Utils.summarizeRun(workflow, params, log)

println params

class ParamsCheck {
    def params
    def definitions

    ParamsCheck(def params, def definitions) {
        this.params = params
        this.definitions = definitions
    }

    def validateAndBuild() {
        def validated = [:]
        def config = this

        definitions.each { paramName, defs ->
            def value = config.params[paramName]

            // Error if required but not provided
            if (defs.required && (value == null || (value instanceof String && value.trim().isEmpty()))) {
                error("Parameter '${paramName}' is required but was not provided.")
            }

            // Handle default value if no value is given
            if (value == null && defs.default_value != null) {
                value = defs.default_value
            }

            // If still no value, continue to next param
            if (value == null) {
                return
            }

            // Trim string values
            if (value instanceof String) {
                value = value.trim()
            }

            // Validate against 'allow' regex if defined
            if (defs.allow) {
                def matched = defs.allow.any { it == '*' || value.matches(it) }
                if (!matched) {
                    error("Value for parameter '${paramName}' does not match allowed patterns: '${value}'.")
                }
            }

            // Convert value to correct type if possible
            if (defs.type == 'integer') {
                try {
                    value = value as int
                } catch (Exception e) {}
            } else if (defs.type == 'float') {
                try {
                    value = value as float
                } catch (Exception e) {}
            } else if (defs.type == 'path') {
                def file = new File(value)
                if (!file.exists()) {
                    error("Path for parameter '${paramName}' does not exist: '${value}'.")
                }
                value = file.toString()
            }
            validated[paramName] = value
        }
        return validated
    }
}

def paramChecker = new ParamsCheck(params, params.definitions)
def validated_params = paramChecker.validateAndBuild()

println validated_params
/*
workflow {
    run_fastp(fastq_ch)
    bwa_align(run_fastp.out.trimmed_ch.combine(ref_ch))
    add_read_groups(bwa_align.out.sorted_bam_ch)
    mark_duplicates(add_read_groups.out.rg_bam_ch)
}
*/

