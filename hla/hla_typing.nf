#!/usr/bin/env nextflow

/***************************************************
*  HLA Typing pipeline
*/

// Parse command line options, set other parameters. 
parseParams()


bamchannel = Channel.fromFilePairs("${params.bamDir}*.{bam,bai}"){file->file.name.replaceAll(/.bam|.bai$/,'')}

println("\nInput channel with: ")
bamchannel.combine([file(params.referenceDir)]).view().set{dv_ch} // add the reference dir to every tuple




process deepvariant{
    cpus 48
    memory '48 GB'
    container "google/deepvariant"

    publishDir "${params.outDir}",
	mode: 'copy'

    input:
    tuple sampleID,path(bamfiles),path(referencedir) from dv_ch

    output: 
    path "*vcf.gz" into dvout_ch

    script:
    bam = bamfiles[0]
    bai = bamfiles[1]

    """ 
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=$referencedir/${params.refName}.fa \
	--reads=$bam \
	--output_vcf=${sampleID}.vcf.gz \
	--output_gvcf=${sampleID}.g.vcf.gz\
	--intermediate_results_dir /intermediate_results_dir \
	--num_shards=48    
    """
}    


process phase_variants{

}


process altref{


}


process calltype{

}


