#!/usr/bin/env nextflow

parseParams()

// It's possible to have any number of files in the "pairs" list.  Default is 2, but with size: can set 
// any size, and with -1 can have unlimited. 
// NOTE: fromFilePairs returns the results in lexicographical order!  
// This kind of sucks... I'd like the pairs returned in the order requested somehow...
Channel.fromFilePairs("${params.inDir}/**.{mcool,pairs.sorted.txt.gz,cool}",size: 3).
set{samples_ch}


process foo{
    input:
    tuple sampleId,file(x) from samples_ch
    
    output:
    stdout out_ch
    
    script:
    println("X: ${x}")
    mcool = x[0]
    pairs = x[1]
    cool =  x[2]
    println("mcool: ${mcool}\t pairs:${pairs}\tcool:${cool}")
    
    """
    echo $x
    """
}
out_ch.view()



def parseParams(){
helpmsg=
"""
    HiNT Nextflow pipeline. 
    
    Run command like:

        hint.nf --inDir ~/coolAndPairsRoot/ --outDir ~/outdir/

    Required arguments:
        --inDir [path]        Name of a directory some subdirectory of which includes .mcool and .pairs.sorted.txt.gz files.
                              NOTE: ALL cool files recursively under this directory will be processed. 
                              
        --outDir [path]       Directory where output files should go.
    
    Optional argument:
        --build               hg19 or hg38 (hg38 default)    
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.inDir || !params.outDir) {log.info helpmsg; exit 1}   

    // Set defaults
    params.build="hg38"
}


