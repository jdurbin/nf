#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process mockbwa {
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
	
    params.outdir = "./tempout/"

    params.inDir = "./testdata/ref/"
    // It's possible to have any number of files in the "pairs" list.  Default is 2, but with size: can set
    // any size, and with -1 can have unlimited.
    // NOTE: fromFilePairs returns the results in lexicographical order!
    // This kind of sucks... I'd like the pairs returned in the order requested. 
    // Also, for some reason -1 returns: 
    //      [bob, [/mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa]]
    //      [bob.fa, [/mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.amb, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.ann, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.bwt, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.pac, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.sa]]
    // 
    // But 5 returns: 
    //      [bob.fa, [/mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.amb, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.ann, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.bwt, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.pac, /mnt/ebs/james/src/nf/sandbox/testdata/ref/bob.fa.sa]]
    //fastq_ch = Channel.fromFilePairs("${params.inDir}/**.{fa,amb,ann,bwt,pac,sa}",size: 5)
    //fastq_ch = Channel.fromFilePairs("${params.inDir}/**.fa*",size: -1)


    // params.fastqDir = "./testdata/fastqs/"
    // fromFilePairs just emits the grouping key with each pair.  Also, a "pair" can be more than two
    //fastq_ch = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
    //            chekIfExists: true,flat: true) 

    fastq_ch | view
    //fastq_ch | mockbwa | view

}
