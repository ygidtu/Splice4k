package com.splice4k.cli


import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.splice4k.index.SJIndex
import com.splice4k.isoforms.Uniquor
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.30
 * @version 20180930
 *
 * 从PacBio中提取构建isoforms
 */


class Iso: CliktCommand(help = "Construct Isoforms through SMRT-seq data") {

    private val input by argument(
            help = "Path to input BAM/SAM files, multiple files separate by space"
    ).file(exists = true).multiple()

        private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file().required()


    private val overlapOfRefReads by option(
            "--overlap-ref-reads",
            help = "Minimal overlap level to match reference with reads [default: 90.0]"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }


    private val show by option(
            "-v",
            "--verbose",
            help = "Enable detailed messages"
    ).flag(default = false)


    override fun run() {

        this.output.apply {
            if ( this.isDirectory ) {
                println("Please set path of output file [event not exists]")
                exitProcess(0)
            }
        }

        val bam = this.input.map {
                SJIndex(
                    infile = it,
                    filter = 0,
                    silent = !this.show,
                    smrt = true
                )
        }

        val uniq = Uniquor(
                reads = bam.first().reads
        )

        uniq.writeTo(this.output)

    }
}