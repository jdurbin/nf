#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process sayHello {
    executor = 'local'

    input:
        val greetings

	output:
		stdout

    script:
        """
        echo '$greetings world!'
        """
}

workflow{
	Channel.from('Bonjour', 'Ciao', 'Hello', 'Hola') \
    | sayHello \
    | map{it.toUpperCase()} \
    | view
}