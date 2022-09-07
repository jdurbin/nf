#!/usr/bin/env nextflow

refdir="s3://dtgb.io/james/hint/ref"
matrixdir="s3://dtgb.io/james/hint/matrix"
genomesizes="s3://dtgb.io/james/scratch/hg38.genome"

def getID(f){
    return(f.name.toString().take(f.name.toString().indexOf('.')))
}

pairsPattern = "/mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/*.pairs.sorted.txt.gz"
coolPattern = "/mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/*.mcool"

pairs = Channel.fromPath(pairsPattern)
    .map{f->id = getID(f); return(tuple(id,f))}


cools = Channel.fromPath(coolPattern)
    .map{f->id = getID(f); return(tuple(id,f))}

pairs.join(cools)
.view() 
.set{samples_ch}

process foo{
    input:
    tuple sampleId,path(pairs),path(mcool) from samples_ch  // note path is newer and prefered to file
    
    output:
    stdout out_ch

    script:
    
    println("pairs: ${pairs} mcool: ${mcool}")   

    """
    echo runsomecommand $pairs $mcool
    cat $genome
    """
}
out_ch.view()
