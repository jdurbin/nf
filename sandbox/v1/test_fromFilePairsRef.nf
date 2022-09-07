#!/usr/bin/env nextflow

parseParams()

println("INDIR: "+params.inDir)

// It's possible to have any number of files in the "pairs" list.  Default is 2, but with size: can set 
// any size, and with -1 can have unlimited. 
// NOTE: fromFilePairs returns the results in lexicographical order!  
// This kind of sucks... I'd like the pairs returned in the order requested somehow...
Channel.fromFilePairs("${params.inDir}/**.{fa,amb,ann,bwt,pac,sa}",size: -1)
.view()
.set{samples_ch}


process foo{
    input:
    tuple sampleId,file(x) from samples_ch
    
    output:
    stdout out_ch
    
    script:
    println("X: ${x}")
    filea = x[0]
    fileb = x[1]
    filec =  x[2]
    println("filea: ${filea}\t fileb: ${fileb}\tfilec: ${filec}")
    
    """
    echo $x
    """
}
out_ch.view()



def parseParams(){
println("PARSE")	
	
helpmsg=
"""
    Required arguments:
        --inDir [path]        Name of a directory some subdirectory of which includes .mcool and .pairs.sorted.txt.gz files.
                              NOTE: ALL cool files recursively under this directory will be processed. 
                              
        --outDir [path]       Directory where output files should go.
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.inDir || !params.outDir) {log.info helpmsg; exit 1}   

println("${params.inDir}")
}


