nextflow.enable.dsl=2

params.trimmed_dir = ""
params.reference = ""
params.out_dir = "."

/*
trimmed_fq_ch = Channel.fromFilePairs("${params.trimmed_dir}/*_{1,2}.fastq", flat: true)
ref_ch = Channel.fromPath("${params.reference}/*")
         .collect()
         .map {
             tuple getCommonPrefix(it), it
         }                                       

def getCommonPrefix(filenames) {
    if (filenames.isEmpty()) return ""
    def prefix = filenames[0].getFileName().toString()
    filenames[1..-1].each { file ->
        file = file.getFileName().toString()
        def i = 0 
        while (i < Math.min(prefix.length(), file.length()) && prefix[i] == file[i]) {
            i++
        }
        prefix = prefix[0..<i]
    }
    return prefix[0..-2] // last common character should always be a period. 
}

input_ch = trimmed_fq_ch.combine(ref_ch)
*/
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
