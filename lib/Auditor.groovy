import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardOpenOption

class Auditor {
    def logFile
    def manifest
    def userName

    /**
     * Constructor for the Auditor class.
     * Initializes the log file and writes the header.
     * @param vparams - The validated parameters map.
     */
    Auditor(def workflow, Map vparams = [:]) {
        this.manifest = workflow.manifest
        this.userName = workflow.userName ?: "unknown"
        def pipelineName = this.manifest.pipeline.name

        def outDir = vparams.containsKey('out_dir') ? vparams.out_dir.toString() : "logs"
        def outDirPath = Paths.get(outDir)
        if (!Files.exists(outDirPath)) {
            Files.createDirectories(outDirPath)
        }

        def timestamp = new Date().format("yyyyMMdd_HHmmss")
        this.logFile = Paths.get(outDirPath.toString(), "${pipelineName}_${timestamp}.log")

        // Write the initial header to the log file, including a readable manifest
        def header = """
Date: ${new Date()}
User: ${this.userName}
---
Manifest:
  Pipeline Name: ${this.manifest.pipeline.name}
  Author: ${this.manifest.author ?: 'N/A'}
  Version: ${this.manifest.version ?: 'N/A'}
  Description: ${this.manifest.description ?: 'N/A'}
"""
        Files.write(this.logFile, header.bytes, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
    }

    /**
     * Adds an ASCII art file to the very top of the log.
     * @param asciiArtFilePath - The path to the ASCII art text file.
     */
    void addAsciiArt(String asciiArtFilePath) {
        if (Files.exists(Paths.get(asciiArtFilePath))) {
            def asciiArt = Files.readAllBytes(Paths.get(asciiArtFilePath))
            def currentLogContent = Files.readAllBytes(this.logFile)
            def newLogContent = new byte[asciiArt.length + currentLogContent.length]
            System.arraycopy(asciiArt, 0, newLogContent, 0, asciiArt.length)
            System.arraycopy(currentLogContent, 0, newLogContent, asciiArt.length, currentLogContent.length)
            Files.write(this.logFile, newLogContent, StandardOpenOption.TRUNCATE_EXISTING)
        } else {
            println "Warning: ASCII art file not found at ${asciiArtFilePath}"
        }
    }

    /**
     * Writes a new section to the log file.
     * @param sectionTitle - The title of the section.
     */
    private void writeSection(String sectionTitle) {
        def separator = "-" * 40
        def content = """
${separator}
${sectionTitle}
${separator}
"""
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
    }

    /**
     * Prints a map of parameters in a human-readable way.
     * @param vparams - The map of parameters.
     */
    void addRunParameters(Map vparams) {
        writeSection("Run Parameters")
        def content = vparams.collect { key, value -> "  ${key}: ${value}" }.join("\n") + "\n"
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
    }

    /**
     * Adds pipeline inputs to the log in a readable format.
     * @param inputs - The inputs to log. This can be a map, list, or string.
     */
    void addPipelineInputs(def inputs) {
        writeSection("Pipeline Inputs")
        def content = inputs.toString() + "\n"
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
    }

    /**
     * Reads a multi-line string and prints it in a formatted "Pipeline Metadata" section.
     * @param metadata - The multi-line string to print.
     */
    void addPipelineMetadata(String metadata) {
        writeSection("Pipeline Metadata")
        def content = metadata + "\n"
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
    }

    /**
     * Appends a new entry to the "Pipeline Outputs" section.
     * @param data - The data to add. This can be a map, list, or string.
     */
    void addPipelineOutput(def data) {
        writeSection("Pipeline Outputs")
        def content = ""
        if (data instanceof Channel) {
            data.subscribe { item ->
                content += "  - ${item}\n"
            }
        } else if (data instanceof List || data instanceof Map) {
            content = data.collect { k, v -> "  ${k}: ${v}" }.join("\n") + "\n"
        } else {
            content = data.toString() + "\n"
        }
        Files.write(this.logFile, content.bytes, StandardOpenOption.APPEND)
    }

    /**
     * Logs the completion status and time.
     */
    void logCompletion() {
        def completionLog = """
----------------------------------------
Pipeline Status: Success
Completion Time: ${new Date()}
----------------------------------------
"""
        Files.write(this.logFile, completionLog.bytes, StandardOpenOption.APPEND)
    }

    /**
     * Logs the failure status and the error that occurred.
     * @param error - The Throwable or error message from the failed pipeline.
     */
    void logFailure(def error) {
        def errorLog = """
----------------------------------------
Pipeline Status: Failure
Failure Time: ${new Date()}
Error: 
${error.toString()}
----------------------------------------
"""
        Files.write(this.logFile, errorLog.bytes, StandardOpenOption.APPEND)
    }
}