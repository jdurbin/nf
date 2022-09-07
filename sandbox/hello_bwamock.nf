#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process bwa {
    executor = 'local'

    input:
        tuple val(id),file(R1s),file(R2s) 

	output:
		stdout

    script:
        """
        echo 'bwa ${R1s} ${R2s}  output: ${id}.bam'
        """
}



workflow{
	params.fastqDir = "./testdata/fastqs/"
    params.outdir = "./tempout/"

    // fromFilePairs just emits the grouping key with each pair.  Also, a "pair" can be more than two
    fastq_ch = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
                chekIfExists: true,flat: true) 

    //fastq_ch | view
    fastq_ch | bwa | view

}
