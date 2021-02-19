#!/usr/bin/env nextflow


// This style of direct path manipulation sort of violates the spirit of nextflow
// The benefit is that I'm not copying all of the files from EFS to S3, from S3 to an EC2 instance, 
// and the results back again to EFS.  Those three copies would probably take as much or more 
// wall clock time as doing the sort.  nextflow is still doing a lot to make this a fairly 
// simple process that can scale to do lots of work.  

def helpMessage(){
log.info"""

Sort and index all of the bam files found in a shared efs path. 

Usage: 

sortandindex.nf --bampath /efs/ref_push/Dev/

--bampath       All bamfiles in this path will be sorted and indexed.  
                WARNING:  This assumes that the path is mounted on batch EC2 instances. 

"""
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
} 

// Use some plain old Groovy to generate a list of the files we want to process...
dir = new File(params.bampath)
fileList=[]
dir.eachFileMatch(~/.*bam/){f-> 
    if (f.name.contains(""))
    fileList.add([f.parent,f.name.replace(".bam","")])
}

// Give the user a heads up on what files are being processed.
println("The following files will be sorted and indexed: ")
fileList.each{
    println(it)
}

// Create a channel from this list...
Channel.fromList(fileList).set{bam_ch}


process sortAndIndexBAM{
    cpus = 8
    memory = '20G'    // I'm requesting total memory for all 8 processes... or should this be per vcpu?
    executor = 'awsbatch'
    queue = 'nfq-VIP'
    container '916851041342.dkr.ecr.us-west-2.amazonaws.com/pyselva' // This works with nfq-VIP queue
    
    input:
    tuple val(bampath),val(bamID) from bam_ch
    
    output:
    stdout out_ch
    """
    samtools sort -m 2G -@ 8 ${bampath}/${bamID}.bam -T /mnt/ebs/tmpjd/${bamID} -o ${bampath}/${bamID}_sorted.bam
    samtools index ${bampath}/${bamID}_sorted.bam
    """
}

out_ch.view()



