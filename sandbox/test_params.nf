#!/usr/bin/env nextflow

parseParams()

println("graph dir: ${params.graphDir}")

cheers = Channel.from 'Bonjour', 'Ciao'

process hla_la_single {
    echo true
    input: 
    val x from cheers
    
    script:
    """
    echo '${x}!  Tachyon beam fired!'
    """
}

/*************************
*   Helper Functions
**/
def parseParams(){
helpmsg=
"""
    foo --phaseVariance 0.00123 
    
    Run command like:
    
        ~/nf/sandbox/test_params.nf 

    Required arguments:
        --phaseVariance         Phase variance for coherent tachyon beam. 
        --graphDir [path]       Graph directory (ideally s3 bucket)   
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.phaseVariance) {log.info helpmsg; exit 1} 
    
    // Default parameters.   
    // It's kind of vague when nextflow replaces parameter assignments 
    // in the script with command line parameters but this works. 
    // This means that it executes parseParams before it subs in the 
    // command line parameter, which is kind of weird. Like, parseParams
    // is good enough to set the default and prevent triggering an error
    // because variable is unset, but somehow it 
    params.graphDir = 's3://dtgb.io/james/scratch/graphs'
}

// This doesn't work, we get a param not defined error in the first
// println.  On the other hand, putting it in parseParams does work. 
// This implies that execution order matters and it's not just 
// prescreening the file for parameter assignments. 
// params.graphDir = 's3://dtgb.io/james/scratch/graphs'
// 
// But if execution order matters, when does the command line verison
// override the default?  
