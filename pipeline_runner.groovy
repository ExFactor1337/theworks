import ParamsChecker
import groovy.util.ConfigSlurper

def banner_art = "assets/banner.txt" // Replace with your actual path
// --- 1. Load Parameter Definitions ---
def configFile = new File("configs/base.config")
def config = new ConfigSlurper().parse(configFile.toURL())
def definitions = config.params.definitions

// --- 2. Convert CLI Arguments to a Map ---
// Process the 'args' array into a map that ParamsChecker expects.
def userParams = [:]
def lastKey = null
args.each { arg ->
    if (arg.startsWith('--')) {
        def parts = arg.substring(2).split('=', 2)
        def key = parts[0]
        
        if (parts.size() > 1) {
            userParams[key] = parts[1]
            lastKey = null
        } else {
            userParams[key] = true
            lastKey = key
        }
    } else if (lastKey != null) {
        userParams[lastKey] = arg
        lastKey = null
    }
}

// --- 3. FAIL-FAST CHECK: Check for help flag and exit successfully ---
def checker = new ParamsChecker(definitions)

// Check for both the key 'help' being present AND it having a truthy value (or just being present)
if (userParams.containsKey('help')) {
    // If the user provided --help, print the guide and exit cleanly.
    checker.printHelp(banner_art) 
    checker.printHelp()
    System.exit(0)
}

// --- 4. Validate Parameters and Get Final Config ---
def vparams

try {
    // If no help flag, proceed with validation of all parameters
    vparams = checker.validate(userParams)
} catch (Exception e) {
    // Print the error, then the help message, then exit with failure code.
    println "\n********************************************************"
    println "ERROR: Parameter validation failed."
    println "${e.getMessage()}"
    println "********************************************************"
    
    // Print the help message for context
    checker.printHelp(banner_art)     
    System.exit(1)
}

// --- 5. Build Nextflow Command ---
def nfCmd = "nextflow run theworks.nf"

// Use the validated parameters (vparams) to build the command string.
vparams.each { k, v -> 
    if (v instanceof Boolean) {
        if (v) {
            nfCmd += " --${k}"
        }
    } else {
        // REMOVE THE SINGLE QUOTES. Let Groovy's execute handle tokenization.
        // If the value is a directory path, passing it unquoted is usually best.
        // If quoting is absolutely needed, use double quotes, and handle internal single quotes.
        
        // Safest approach: pass arguments without literal single quotes unless truly necessary.
        nfCmd += " --${k} ${v}" 
    }
}

println "Running command:\n$nfCmd"


// Execute the command and handle errors
def proc = nfCmd.execute()

proc.in.eachLine { println it }
proc.err.eachLine { System.err.println it }
proc.waitFor()
if (proc.exitValue() != 0) {
    println "Nextflow run failed with exit code ${proc.exitValue()}"
    System.exit(proc.exitValue())
} else {
    println "Nextflow run completed successfully."
}
