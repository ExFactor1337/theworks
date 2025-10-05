import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardOpenOption

class Auditor {
    def logFile
    def manifest
    def userName

    // ... (Constructor and other methods remain the same) ...
    Auditor(def workflow, Map vparams = [:]) {
        this.manifest = workflow.manifest
        this.userName = workflow.userName ?: "unknown"
        def pipelineName = this.manifest.name

        def outDir = vparams.containsKey('out_dir') ? vparams.out_dir.toString() : "logs"
        def outDirPath = Paths.get(outDir)
        if (!Files.exists(outDirPath)) {
            Files.createDirectories(outDirPath)
        }

        def timestamp = new Date().format("yyyyMMdd_HHmmss")
        this.logFile = Paths.get(outDirPath.toString(), "${pipelineName}_${timestamp}.log")
    }

    // Existing instance method for logging (modified to delegate print to static method)
    void addAsciiArt(String asciiArtFilePath) {
        if (Files.exists(Paths.get(asciiArtFilePath))) {
            def asciiArt = Files.readAllBytes(Paths.get(asciiArtFilePath))
            // Write to log file
            Files.write(this.logFile, asciiArt, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
            // Print to stdout
            println new String(asciiArt) 
        } else {
            println "Warning: ASCII art file not found at ${asciiArtFilePath}"
        }
    }

    // NEW STATIC METHOD: Prints the ASCII banner without needing a log file or full setup
    static void printAsciiBanner(String asciiArtFilePath) {
        if (Files.exists(Paths.get(asciiArtFilePath))) {
            def asciiArt = Files.readAllBytes(Paths.get(asciiArtFilePath))
            println new String(asciiArt)
        } else {
            println "Warning: ASCII art file not found at ${asciiArtFilePath}"
        }
    }

    // ... (Other methods: printBasicInfo, printVparams, etc. remain the same) ...
    void printBasicInfo(Map vparams) {
        def info = """
========================================
Pipeline Name: ${this.manifest.pipeline.name} ?: 'Polly Pipeline'
Nextflow Scratch Dir: ${System.getenv('NXF_WORK')}
Run Directory: ${System.getProperty('user.dir')}
User: ${this.userName}
Start Time: ${new Date()}
========================================
"""
        Files.write(this.logFile, info.bytes, StandardOpenOption.APPEND)
        println info

        // Help info
        def help = "For help, see documentation or contact pipeline maintainer."
        Files.write(this.logFile, ("\n" + help + "\n").bytes, StandardOpenOption.APPEND)
        println help
    }

    void printVparams(Map vparams) {
        def content = "\nParameters:\n" + vparams.collect { key, value -> "  ${key}: ${value}" }.join("\n") + "\n"
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
        println content
    }

    void printOutputFiles(List outputFiles) {
        def content = "\nOutput Files:\n" + outputFiles.collect { f -> "  - ${f}" }.join("\n") + "\n"
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
        println content
    }

    void printTimestamp() {
        def ts = "\nTimestamp: ${new Date()}\n"
        Files.write(this.logFile, ts.bytes, StandardOpenOption.APPEND)
        println ts
    }

    void printStatus(String status, def error = null) {
        def statusLine = "\nStatus: ${status}\n"
        Files.write(this.logFile, statusLine.bytes, StandardOpenOption.APPEND)
        println statusLine
        if (status == "Failure" && error != null) {
            def errorMsg = "Error:\n${error.toString()}\n"
            Files.write(this.logFile, errorMsg.bytes, StandardOpenOption.APPEND)
            println errorMsg
        }
    }
    void generateAudit(Map vparams, List outputFiles, String asciiArtFilePath) {
        this.addAsciiArt(asciiArtFilePath)
        this.printBasicInfo(vparams)
        this.printVparams(vparams)
        this.printOutputFiles(outputFiles)
        this.printTimestamp()
    }
}