#!/usr/bin/env nextflow


def parseParams(){
helpmsg=
"""
    Cantata Bio HLA Typing and Phasing
    
    Run command like:
    
    hla_typephase.nf --fastqdir myfastqs --outDir myoutdir

    Required arguments:
        --fastqDir [path]       Name of fastq dir.         
        --outDir [path]         Directory where output files should go.         
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.fastqDir || !params.outDir) {        
        log.info helpmsg; exit 1
    }     
}
parseParams()

Channel
.fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",chekIfExists: true,flat: true)
.set{fastqDir_ch}	

process bwa_mock{
    input:
    tuple sampleId,file(x),file(y) from fastqDir_ch

    output: 
        stdout out_ch 

    script:
        println("X: ${x} Y: ${y}") 

        """
        echo runsomecommand $x $y
        """

}

out_ch.view()