#!/usr/bin/env nextflow

import grapnel.util.Parser

//println("Args are: ${args}")

println(params)


pargs = params.collect{key,value->[key,value]}.flatten()
println(pargs)

def ParseOptions(params){
    pargs = params.collect{key,value->[key,value]}.flatten()
    
	parser = new Parser(description: '''
    Nextflow pipeline for running HiNT    
    ''');   
    try{
        parser.with{
            required 'c','coolDir',[description:'Directory with bam files']
            required 'o','outDir',[description: 'Directory where output files should go.']
            optional 'b','build',[default:'hg38',description: 'Optional genome build.',validate:{(it=='hg38' || it=='hg19')}]
        }
        options = parser.parse(pargs); 
    }catch(Exception e){
        System.err.println "Exception :"+e
        System.err << parser.usage;
        System.exit(1)
    }		
	return(options)
}

params = ParseOptions(params)