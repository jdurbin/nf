#!/usr/bin/env nextflow

params.mcoolDir = false
params.outDir = false
params.resolutions = false
params.help   = false

#!/usr/bin/env nextflow

params.mcoolDir = false
params.outDir = false
params.resolutions = false
params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        hic_breakfinder.nf --bamDir ~/mybamfiledir/

    Mandatory arguments:
        --bamDir [path]           Name of the a directory with bam files. 
        --outDir [path]           Directory where output files should go.
    
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.bamDir) {
        exit 1, "--bamDir is a required argument. Use --help to get the full usage." 
}

Channel
    .fromPath("${params.bamDir}/*.bam")
    .set{hicbreakfinder_ch}

process hic_breakfinder {
    cpus 8
    memory '48 GB'
    container "kjdurbin/hic_breakfinder"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    path(bam) from hicbreakfinder_ch
    
    output:
    tuple id,path("*.txt") into hicbreakfinder_out_ch
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    /usr/local/bin/hic_breakfinder --bam-file ${bam} \
    --exp-file-inter /src/hic_breakfinder/resources/inter_expect_1Mb.hg38.txt \
    --exp-file-intra /src/hic_breakfinder/resources/intra_expect_100kb.hg38.txt --name ${id}    
    """
}