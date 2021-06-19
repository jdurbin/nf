#!/usr/bin/env nextflow

/***************************************************
*  HiNT translocation finding pipeline
*/

// Parse command line options, set other parameters. 
parseParams()

// Create channels from inputs
bamchannel = Channel.fromPath(params.bamFile).map{f->id = getID(f); return(tuple(id,f))}

println("\nInput channel with: ")
bamchannel.combine([file(params.referenceDir)]).view().set{dv_ch} // add a singleton

// Process
process dv_simple {
    cpus 48
    memory '48 GB'
    container "google/deepvariant"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple sampleID,path(bamfile),path(referencedir) from dv_ch
    
    output: 
    path "*vcf.gz" into dvout_ch
    
    script:
    
    """ 
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=$referencedir/hg38.fa \
	--reads=$bamfile \
	--output_vcf=${sampleID}.vcf.gz \
	--output_gvcf=${sampleID}.g.vcf.gz\
	--intermediate_results_dir /intermediate_results_dir \
	--num_shards=48    
    """
}


/*************************
*   Helper Functions
**/

def parseParams(){
helpmsg=
"""
    Simple deepvariant pipeline
    
    Run command like:
    
        ~/nf/deepvariant/dv_simple.nf --bam myfile.bam --outDir dvout/ --referenceFile hs37d5.fa.gz

    Required arguments:
        --bamFile                 Bamfile for variant calling.  
        --outDir [path]           Directory where output files should go.
        --referenceDir [path]     Directory containing reference genome and associated indices. 
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.bamFile || !params.outDir || !params.referenceDir) {log.info helpmsg; exit 1} 
}


// Extract the portion of the name up to the first dot to use as an ID
def getID(f){
    return(f.name.toString().take(f.name.toString().indexOf('.')))
}
