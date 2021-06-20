#!/usr/bin/env nextflow


// It's possible to have any number of files in the "pairs" list.  Default is 2, but with size: can set 
// any size, and with -1 can have unlimited. 
// NOTE: fromFilePairs returns the results in lexicographical order!  
// This kind of sucks... I'd like the pairs returned in the order requested somehow...
//Channel.fromFilePairs("/mnt/ebs/james/src/nf/sandbox/testdata/bamfiles/*").view().set{samples_ch}

//Channel
//    .fromFilePairs('/mnt/ebs/james/src/nf/sandbox/testdata/bamfiles/*.{bam,bai}') { file -> file.name.replaceAll(/.bam|.bai$/,'') }
//    .view()
//    .set { samples_ch }


Channel
    .fromFilePairs('/mnt/ebs/james/src/nf/sandbox/testdata/bamfiles/*.{bam,bai}'){file->file.name.replaceAll(/.bam|.bai$/,'')}
    .view()
    .set { samples_ch }



//Channel.fromPath("/Users/james/src/nf/sandbox/testdata/bamfiles/*").view().set{samples_ch}


process foo{
    input:
    tuple sampleID,file(bslist) from samples_ch
    
    output:
    stdout out_ch
    
    script:
    bam = bslist[0]
    bai = bslist[1]
   
    """
    echo bam: ${bam}
    echo bai: ${bai}
    
    cat $bam
    cat $bai
    """
}
out_ch.view()