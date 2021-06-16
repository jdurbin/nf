#!/usr/bin/env nextflow

params.coolDir = false
params.outDir = false
params.help   = false
params.build="hg38"

parseParams()

Channel.fromPath('/mnt/ebs/james/dan/hint/nf/testdata/**.mcool').set{samples_ch}

process foo{
    input:
    file x from samples_ch
    
    output:
    stdout out_ch
    
    script:
    """
    echo $x
    """
}
out_ch.view()




/**
*   Parse options and give help
*/ 
def parseParams(){
helpmsg=
"""
    HiNT Nextflow pipeline. 
    
    Run command like:

        hint.nf --coolDir ~/mycool/ --outDir ~/outdir/

    Required arguments:
        --coolDir [path]        Name of the a directory with bam files. 
        --pairsDir [path]       Name of the a directory with pairs files. 
        --outDir [path]         Directory where output files should go.
    
    Optional argument:
        --build                   hg19 or hg38 (hg38 default)    
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    //if (!params.coolDir || !params.outDir || !params.pairDir) {log.info helpmsg; exit 1}    
}