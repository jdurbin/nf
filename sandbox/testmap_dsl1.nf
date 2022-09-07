#!/usr/bin/env nextflow
//nextflow.enable.dsl=2

params.fastqDir="./testdata/fastqs/"

Channel
.empty()
.into {fastq_ch}

fastqDir_ch = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
                chekIfExists: true,flat: true) 


process mockbwa {
    executor = 'local'

    input:
        tuple id, file(R1s), file(R2s) from fastq_ch
        .map{bs,file ->
	        pref = file.name.toString().take(file.name.toString().indexOf('_R'))
	        return(tuple(pref,file))
	    }
	    .groupTuple()
	    .flatten()
	    .collate(3)
        .mix(fastqDir_ch)

	output:
		stdout outch

    script:
        """
        echo 'bwa ${R1s} ${R2s}  output: ${id}.bam '
        """
}

outch.view()

