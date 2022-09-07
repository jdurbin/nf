#!/usr/bin/env nextflow

cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

process sayHello {
    executor = 'local'

    input:
        val x from cheers

	output:
		stdout out_ch

    script:
        """
        echo '$x world!'
        """
}

out_ch.view()