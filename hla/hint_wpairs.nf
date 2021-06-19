#!/usr/bin/env nextflow

/***************************************************
*  HiNT translocation finding pipeline
*/

// Parse command line options, set other parameters. 
parseParams()

// Define inputs
pairsPattern = "${params.pairsDir}*.pairs.sorted.txt.gz"
coolPattern = "${params.coolDir}/*.mcool"

// Create channels from inputs
pairs = Channel.fromPath(pairsPattern).map{f->id = getID(f); return(tuple(id,f))}
cools = Channel.fromPath(coolPattern).map{f->id = getID(f); return(tuple(id,f))}
joinchannel = pairs.join(cools) // joins channels based on id

println("\nInput channel with: ")
joinchannel.combine([file(params.supportfiles)]).view().set{hint_ch} // add a singleton

// Process
process hint {
    cpus 1
    memory '48 GB'
    container "kjdurbin/hint"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple sampleId,path(pairsfile),path(mcool),path(supportfiles) from hint_ch
    
    output: 
    path "hintout_*" into hintout_ch
    
    script:
    id = mcool.name.toString().take(mcool.name.toString().lastIndexOf('.'))
    
    """ 
    hint tl -m \
    ${mcool}::/resolutions/1000000,${mcool}::/resolutions/100000 \
    -f cooler \
    --chimeric $pairsfile \
    --refdir ${supportfiles}/ref/hg38 \
    --backdir ${supportfiles}/matrix/hg38 \
    -g hg38 \
    -n ${id} \
    -c 0.05 \
    --ppath /pairix/bin/pairix -p 12 \
    -e DpnII \
    -o hintout_${id}
    """
}


/*************************
*   Helper Functions
**/

def parseParams(){
helpmsg=
"""
    HiNT Nextflow pipeline. 
    
    Run command like:
    
        ~/nf/hint/hint_wpairs.nf -profile local --coolDir cools/ --pairsDir pairs/ --outDir hintout/

    Required arguments:
        --coolDir [path]        Name of the a directory with bam files. 
        --pairsDir [path]       Name of the a directory with pairs files. 
        --outDir [path]         Directory where output files should go.
    
    Optional argument:
        --build                 hg19 or hg38 (hg38 default)
        --supportfiles          s3 path to hint support files (ref and matrix directories)    
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.coolDir || !params.outDir || !params.pairsDir) {log.info helpmsg; exit 1} 
    params.build="hg38"
    params.supportfiles="s3://dtgb.io/james/hint/"
}


// Extract the portion of the name up to the first dot to use as an ID
def getID(f){
    return(f.name.toString().take(f.name.toString().indexOf('.')))
}
