import java.nio.file.Paths
import java.text.SimpleDateFormat
class Utils {
    public static String getCommonPrefix(filenames){
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
        return prefix[0..-2]
    }
    public static void summarizeRun(workflow = workflow, params = params, log = log){
        def summary = [:]
        summary['Pipeline Name'] = 'theWorks'
        summary['Version'] = '1.0'
        summary['Config Files'] = workflow.configFiles.join(', ')
        summary['User'] = workflow.userName
        summary['Start Time'] = this.getTimestamp()
        if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - ENABLED"
        summary['Launch Dir'] = workflow.launchDir
        summary['Project Dir'] = workflow.projectDir
        summary['Working Dir'] = workflow.workDir
        summary['Output Dir'] = params.out_dir
        log.info "--------------------SUMMARY--------------------"
        log.info summary.collect {
            k,v -> if(k != null || v != null) "${k.padRight(18)}: $v" 
            }.join("\n") 
        log.info "-----------------------------------------------"
    }
    public static String getTimestamp(){
        def timestamp = new SimpleDateFormat("yyyy-MM-dd")
        return timestamp.format(new Date())
    }
    public static String toAbsPath(stringPath){ 
        def path = Paths.get(stringPath)
        return path.toAbsolutePath().toString()
    }
}
