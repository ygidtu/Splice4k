package com.splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.index.SJIndex

/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180928
 */


class Extract: CliktCommand(help = "Extract splice junctions from BAM/SAM files") {
    private val input by option(
            "-i",
            "--input",
            help = "Path to input BAM/SAM file or gmap align (-A parameter) output file"
    ).file(exists = true).required()

    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file().required()

    private val junctionsFilter by option(
            "-c",
            "--count",
            help = "Filter low abundance junctions [default: 0]"
    ).int().default(0).validate {
        it >= 0
    }


    override fun run() {
        SJIndex(
                infile = this.input.absoluteFile,
                filter = this.junctionsFilter,
                silent = true,
                smrt = false
        ).writeTo(output = this.output.absoluteFile)
    }
}