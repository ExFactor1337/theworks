import java.nio.file.Paths
import java.text.SimpleDateFormat
class Utils {
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
        summary['Output Dir'] = params.outDir
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
