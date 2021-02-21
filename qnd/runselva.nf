#!/usr/bin/env nextflow


def helpMessage(){
log.info"""

Run the selva pipeline on a set of .mcool and .bam files. 

WARNING:  This script assumes that the paths are mounted on 
          all batch EC2 instances. 

Usage: 

runselva.nf \
--coolpath /mnt/ebs/coolerFiles/ \
--bampath /mnt/ebs/sorted_indexed/ \
--outdir /mnt/ebs/selva/output2/RAW/

--coolpath          Selva will be run on all mcool files in this path. 
--bampath           Each .mcool is expected to have a corresponding .bam file with same name and .bam extension. 
--outdir            Directory where output should go. 
--reference         Full path to reference genome (hg38 is default)

"""
}

// Show help message
if (!params.coolpath){
    helpMessage()
    exit 0
} 

// Default genome
params.genome="/mnt/ebs/genome/hg38/hg38.fa"

outdir = new File(params.outdir)
if (!outdir.exists()) throw new Exception("Dir doesn't exist: ${outdir}")


// Use some plain old Groovy to generate a list of the files we want to process...
dir = new File(params.coolpath)
fileList=[]
//dir.eachFileMatch(~/.*mcool/){f-> 
dir.eachFile(){f->
    id = f.name
    coolin = "${f.parent}/${f.name}/${f.name}.mcool"
    // Every .mcool file should have a corresponding sorted bam file...
    bamin = "${params.bampath}/${id}_sorted.bam"
    bamfile = new File(bamin)
    if (!bamfile.exists()){
        throw new Exception("Bamfile doesn't exist: ${bamin}")
    }
    fileList.add([id,bamin,params.outdir,coolin,params.genome])
}

// Give the user a heads up on what files are being processed.
println("Selva will be run on the following contact matrices: ")
fileList[2..-1].each{
    println(it[3])
}


// Create a channel from this list...

Channel.fromList(fileList[2..-1]).set{work_list}

process runSelva {
    cpus = 25
    memory = '100G'   
    executor = 'awsbatch'
    queue = 'nfq-VIP'
    container '916851041342.dkr.ecr.us-west-2.amazonaws.com/pyselva' // This works with nfq-VIP queue

    input:
    tuple id,bamin,outdir,coolin,genome from work_list

    output:
    stdout out_ch
    """
    python3 /selva_tad/bin/cooler_selva_pipe.py find-sv ${id} ${bamin} ${outdir} ${coolin} ${genome} --threads=25
    """
}
out_ch.view()








