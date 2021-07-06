#!/usr/bin/env nextflow

/***************************************************
*  HLA*LA Single BAM pipeline
*/

// Parse command line options, set other parameters. 
parseParams()

// Pick up the bam and associated bai files...
bamRoot = params.bamFile.replace(".bam","")
Channel.fromFilePairs("${bamRoot}*.{bam,bai}"){ 
    file -> file.name.replaceAll(/.bam|.bai$/,'') }.view().set{bamin_ch}

// Add in the pre-computed graph directory (default on s3)
bamin_ch.combine([file(params.graphDir)]).view().set{hlain_ch} // add a singleton


/***************************************************
*  hla_la_single
*/
process hla_la_single {
    echo true
    
    cpus 48
    memory '100 GB'
    container "kjdurbin/hla-la"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple id,path(bamfiles),path(graphdir) from hlain_ch
        
    output: 
    path "working/*" into hlalaout_ch
    stdout out_ch
    
    script:
    bamfile=bamfiles[0]
    baifile=bamfiles[1]
    
    sid = id.replaceAll("-","_") // HLA-LA doesn't like dashes in ID
    
    """
    mkdir -p working
    HLA-LA.pl --bam ${bamfile} --sampleID $sid --workingDir working --customGraphDir ${graphdir} \
    --graph PRG_MHC_GRCh38_withIMGT --maxThreads 48 --picard_sam2fastq_bin '/picard.jar'
    """ 
}
out_ch.view()


/*************************
*   Helper Functions
**/
def parseParams(){
helpmsg=
"""
    HLA*LA Simple Pipeline
    
    Run command like:
    
        ~/nf/hla/hla_la_single.nf -profile local --bamFile input.bam --outDir hintout/

    Required arguments:
        --bamFile [path]        Name of the bam file. 
        --outDir [path]         Directory where output files should go. 
        --graphDir [path]       Graph directory (default: s3://dtgb.io/james/hla-la/graphs/)
""".stripIndent()    
    
    if (params.help){log.info helpmsg; exit 0}
    if (!params.bamFile || !params.outDir) {log.info helpmsg; exit 1} 
    
    // Default graph dir if one isn't given on the command line. 
    // Note: It's kind of weird that this works here, as when does the command line 
    // version override it?   Nonetheless, this does work. 
    params.graphDir = 's3://dtgb.io/james/hla-la/graphs/'
}

// Extract the portion of the name up to the first dot to use as an ID
def getID(f){
    return(f.name.toString().take(f.name.toString().indexOf('.')))
}
