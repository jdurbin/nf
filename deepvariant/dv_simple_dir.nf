#!/usr/bin/env nextflow

/***************************************************
*  HiNT translocation finding pipeline
*/

params.refName=false

// Parse command line options, set other parameters. 
parseParams()

bamchannel = Channel.fromFilePairs("${params.bamDir}*.{bam,bai}"){file->file.name.replaceAll(/.bam|.bai$/,'')}

println("\nInput channel with: ")
bamchannel.combine([file(params.referenceDir)]).view().set{dv_ch} // add a singleton

// Process
process dv_simple_dir {
    cpus 48
    memory '48 GB'
    container "google/deepvariant"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple sampleID,path(bamfiles),path(referencedir) from dv_ch
    
    output: 
    path "*vcf.gz" into dvout_ch
    
    script:
    bam = bamfiles[0]
    bai = bamfiles[1]
    
    """ 
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=$referencedir/${params.refName}.fa \
	--reads=$bam \
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
    
        ~/nf/deepvariant/dv_simple.nf --bam mybamdir --outDir dvout/ --referenceFile hs37d5.fa.gz

    Required arguments:
        --bamDir                  Directory containing bam files and associated bai indices.   
        --outDir [path]           Directory where output files should go.
        --referenceDir [path]     Directory containing reference genome and associated indices. 
        --refName                 Root name of reference
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.refName || !params.bamDir || !params.outDir || !params.referenceDir) {log.info helpmsg; exit 1} 
}

