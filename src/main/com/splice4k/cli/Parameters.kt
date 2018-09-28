package com.splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180928
 */


class Parameters: CliktCommand(invokeWithoutSubcommand = true) {
    private val version by option("-v", "--version", help = "version").flag(default = false)

    override fun run() {
        if ( this.version ) {
            println("Splice4k version: 20180928")
            exitProcess(0)
        }
    }
}