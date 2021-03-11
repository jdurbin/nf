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

--pairspath     All pairsfiles in this path will be sorted and indexed.  
                WARNING:  This assumes that the path is mounted on batch EC2 instances. 

"""
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
} 

pairIDs=["CT432","CT433","CT434","CT435","MB290","MSB549","MSB550","MSB551","MSB552","MSB553","MSB554","MSB700","MSB701","MSB703","MSB704"]
//pairIDs=["CT432"]
rootPath="/mnt/ebs/ref_push/Dev/fromSanjit/coolerFiles"

fileList = pairIDs.collect{pairID->
    [rootPath,"${rootPath}/${pairID}/${pairID}.pairs.sorted.txt.gz",pairID]
}

// Give the user a heads up on what files are being processed.
println("The following files will made into .hic files: ")
fileList.each{
    println(it)
}

// Create a channel from this list...
Channel.fromList(fileList).set{pairs_ch}

process makeHiC{
    cpus = 8
    memory = '20G'    // I'm requesting total memory for all 8 processes... or should this be per vcpu?
    //container '916851041342.dkr.ecr.us-west-2.amazonaws.com/pyselva' // This works with nfq-VIP queue
    container 'kjdurbin/pairtools' // This works with nfq-VIP queue
    
    input:
    tuple val(rootPath),val(pairsFile),val(pairsID) from pairs_ch
    cpus = 8
    
    output:
    stdout out_ch
    """
    cut -f 1,2 /mnt/ebs/genome/hg38/hg38.fa.fai > hg38.genome
    chrlength=hg38.genome

    pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' ${pairsFile} -o ${rootPath}/${pairsID}/${pairsID}.filtered.pairs.gz

    zcat ${rootPath}/${pairsID}/${pairsID}.filtered.pairs.gz | grep -v "#" | awk '{print "1\t"\$2"\t"\$3"\t1\t0\t"\$4"\t"\$5"\t0"}' | awk '{if (\$2 > \$6) {print \$1"\t"\$6"\t"\$7"\t"\$8"\t"\$5"\t"\$2"\t"\$3"\t"\$4} else {print}}' | sort -k2,2d -k6,6d  -T ./ --parallel=8 > ${rootPath}/${pairsID}/${pairsID}_juicer_alignments.txt
    
    java -jar -Xmx30g ~/ebs/james/src/dovetail_tools/juicer_tools_1.22.01.jar pre ${rootPath}/${pairsID}/${pairsID}_juicer_alignments.txt ${rootPath}/${pairsID}/${prefix}.hic ${chrlength}
    """
}

out_ch.view()



