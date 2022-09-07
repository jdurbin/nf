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
	def cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

    main: 
	    sayHello(cheers)
	    sayHello.out.view()

    emit: 
        sayHello.out
}