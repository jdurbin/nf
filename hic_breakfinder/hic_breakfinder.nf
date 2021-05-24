#!/usr/bin/env nextflow

params.bamDir = false
params.outDir = false
params.help   = false
params.build="hg38"

def helpMessage() {
    log.info"""
    Run hic_breakfinder on a directory of alignmnt bam files. 

    Run command like:

        hic_breakfinder.nf --bamDir ~/mybamfiledir/ --outDir ~/outdir/

    Mandatory arguments:
        --bamDir [path]           Name of the a directory with bam files. 
        --outDir [path]           Directory where output files should go.
    
    Optional argument:
        --build                   hg19 or hg38 (hg38 default)
    
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.bamDir) {
        exit 1, "--bamDir is a required argument. Use --help to get the full usage." 
}

if (!params.outDir) {
        exit 1, "--outDir is a required argument. Use --help to get the full usage." 
}

if (params.build == "hg19"){
    interExpect = "s3://dtgb.io/james/hic_breakfinder/inter_expect_1Mb.hg19.txt"
    intraExpect = "s3://dtgb.io/james/hic_breakfinder/intra_expect_100kb.hg19.txt"
}else{
    interExpect = "s3://dtgb.io/james/hic_breakfinder/inter_expect_1Mb.hg38.txt"
    intraExpect = "s3://dtgb.io/james/hic_breakfinder/intra_expect_100kb.hg38.txt"
}

// Build a list of tuples for each job...
dir = new File(params.bamDir)
fileList=[]
dir.eachFileMatch(~/.*bam/){f-> 
    id = f.name 
    fileList.add([id,f.path,interExpect,intraExpect])
}

println("Running on the following tuples:")
fileList.each{
    println(it)
}

Channel
.fromList(fileList)
.set{hicbreakfinder_ch}

process hic_breakfinder {
    cpus 8
    memory '48 GB'
    container "kjdurbin/hic_breakfinder"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple val(id),path(bam),path(inter),path(intra) from hicbreakfinder_ch
    
    output:
    tuple id,path("*.txt") into hicbreakfinder_out_ch
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    aws s3 cp ${interExpect} . 
    aws s3 cp ${intraExpect} . 
    
    /usr/local/bin/hic_breakfinder --bam-file ${bam} \
    --exp-file-inter inter_expect_1Mb.hg38.txt \
    --exp-file-intra intra_expect_100kb.hg38.txt --name ${id}
    """
}

// NOTE:  I couldn't get the s3 expect files to automatically resolve like I thought 
// it should.  To get the run done I just manually put in aws cp to copy the expect files
// but I plan to go back and figure out how to do it the nf way. 