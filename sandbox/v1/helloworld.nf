#!/usr/bin/env nextflow

cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

process sayHello {
    executor = 'local'  
    
    input:
        val x from cheers
    
    script:
        """
        echo '$x world!'
        """
}