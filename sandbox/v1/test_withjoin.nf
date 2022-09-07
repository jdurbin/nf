#!/usr/bin/env nextflow

def getID(f){
    return(f.name.toString().take(f.name.toString().indexOf('.')))
}


pairsPattern = "/mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/*.pairs.sorted.txt.gz"
coolPattern = "/mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/*.mcool"

pairs = Channel.fromPath(pairsPattern)
    .map{f->id = getID(f); return(tuple(id,f))}


cools = Channel.fromPath(coolPattern)
    .map{f->id = getID(f); return(tuple(id,f))}
//    .map{f->id = getID(f); return(tuple(id,f))}.view()

/*
    [CT432, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT432.pairs.sorted.txt.gz, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT432.mcool]
    [CT434, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT434.pairs.sorted.txt.gz, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT434.mcool]
    [CT433, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT433.pairs.sorted.txt.gz, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT433.mcool]
    [CT435, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT435.pairs.sorted.txt.gz, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT435.mcool]
*/
pairs.join(cools)
.view() // I like to get a preview of what we are going to get
.set{samples_ch}

/*
[CT432, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT432.mcool, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT432.pairs.sorted.txt.gz]
[CT433, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT433.mcool, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT433.pairs.sorted.txt.gz]
[CT435, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT435.mcool, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT435.pairs.sorted.txt.gz]
[CT434, /mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/CT434.mcool, /mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/CT434.pairs.sorted.txt.gz]
*/
//cools.join(pairs).view().set{samples_ch}


process foo{
    input:
    tuple sampleId,path(pairs),path(mcool) from samples_ch  // note path is newer and prefered to file
    
    output:
    stdout out_ch

    script:
    
    println("pairs: ${pairs} mcool: ${mcool}")   

    """
    echo runsomecommand $pairs $mcool
    """
}
out_ch.view()
