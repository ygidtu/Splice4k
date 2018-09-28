package com.splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180927
 */


class Parameters: CliktCommand(invokeWithoutSubcommand = true) {
    private val version by option("-v", "--version", help = "version").flag(default = false)

    override fun run() { }
}