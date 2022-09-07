#!/usr/bin/env nextflow

//Channel
//    .fromFilePairs("/mnt/ebs/james/src/nf/sandbox/testdata/**.{mcool,pairs.sorted.txt.gz,cool}", size: -1){it->
//        println("it: $it")
//    }
////    .view { ext, files -> "Files with the extension $ext are $files" }
//    .set{samples_ch}
//    
    
//Channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )

//Channel.fromFilePairs("${params.inDir}/**.{mcool,pairs.sorted.txt.gz,cool}",size: 3).
//set{samples_ch}

//Channel
//    .fromFilePairs(
//        ["/mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/*.pairs.sorted.txt.gz",        
//         "/mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/*.mcool"],
//         flat:true    // Default is false and files are just space separated.
//    )
//    .view(){"Channel contains: $it"}
//    .set{samples_ch}

// This works, so far as it goes, but seems to convert file1 and file2 into mere strings
// and so I think it probably isn't marshalling the files
/*
Channel
    .fromFilePairs(
        ["/mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/*.pairs.sorted.txt.gz",        
         "/mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/*.mcool"],
         flat:true  
    ).map{prefix,file1,file2-> tuple(prefix,["file1":file2,"file2":file1])}   
    .view(){"Channel contains: $it"}
    .set{samples_ch}
*/

filePat= ['pairs':"/mnt/ebs/james/src/nf/sandbox/testdata/pairsfiles/*.pairs.sorted.txt.gz",        
          'mcool':"/mnt/ebs/james/src/nf/sandbox/testdata/coolfiles/*.mcool"]

Channel
    .fromFilePairs(
        [filePat['pairs'],filePat['mcool']],
         flat:true    // Default is false and files are just space separated.
    ).map{prefix,file1,file2-> tuple(prefix,file2,file1)} // so this flips them around... 
    .view(){"Channel contains: $it"}
    .set{samples_ch}
    
process foo{
    input:
//    tuple sampleId,val(fileMap) from samples_ch // with flat:false would only get x
    tuple sampleId,file(x),file(y) from samples_ch // with flat:false would only get x as BlankSeparatedList

    output:
    stdout out_ch

    script:
    
    println("X: ${x} Y: ${y}")   

    """
    echo runsomecommand $x $y
    """
}
out_ch.view()