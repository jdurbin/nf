#!/usr/bin/env -S nextflow -log ./work/log/nf.log run 
nextflow.enable.dsl=2

process bwa {
    executor = 'local'

    input:
        tuple val(id),file(R1),file(R2) 

	output:
        tuple(val(id),path("${id}.bam"))

    script:
        """
        echo $R1 $R2 > ${id}.bam
        """
}


process sort{
    executor = 'local'

    input: 
        tuple val(id),file(unsortedbam)

    output:
        path "${id}.sorted.bam"

    script:        
        """
        echo $unsortedbam > ${id}.sorted.bam
        """
}

process mergebam{
    executor = 'local'

    input:
        val(sortedbamfiles)

    output:
        path "merged.bam"

    script:
    // nextflow.util.ArrayBag

        sortedbamStr = sortedbamfiles.join(" ")
        """
        echo $sortedbamStr > merged.bam
        """    
}


workflow{
	params.fastqDir = "./testdata/fastqs/"
    params.bamoutID = "HG002"
    params.outdir = "./tempout/"

    // fromFilePairs just emits the grouping key with each pair.  Also, a "pair" can be more than two
    fastq_ch = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
                chekIfExists: true,flat: true) 

    //fastq_ch | view
    fastq_ch | bwa | sort | collect | mergebam | view

}
