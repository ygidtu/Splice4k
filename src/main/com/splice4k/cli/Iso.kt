package com.splice4k.cli


import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.splice4k.index.SJIndex
import com.splice4k.index.AnnotationIndex
import com.splice4k.isoforms.tools.GeneReadsCoupler
import org.apache.log4j.Logger
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.30
 * @version 20180930
 *
 * 从PacBio中提取构建isoforms
 */


class Iso: CliktCommand(help = "Construct Isoforms through SMRT data") {

    private val input by option(
            "-i",
            "--input",
            help = "Path to input BAM/SAM files, multiple files separate by space"
    ).file(exists = true).required()


    private val reference by option(
            "-r",
            "--reference",
            help = "Path to reference file [gtf|gff3]"
    ).file(exists = true).required()


    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file().required()


//    private val threads by option(
//            "-p",
//            "--process",
//            help = "Number of processes to use"
//    ).int().default(1).validate {
//        require( it > 0 && it <= Runtime.getRuntime().availableProcessors())
//    }
//
//
//    private val error by option(
//            "-e",
//            help = "The error to identify whether AS event exists [default: 3bp]"
//    ).int().default(3).validate {
//        require( it >= 0 ) {"this value must be positive"}
//    }


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
        val logger = Logger.getLogger(Iso::class.java)

        this.output.apply {
            if ( this.isDirectory ) {
                println("Please set path of output file [event not exists]")
                exitProcess(0)
            }
        }


        val ref = AnnotationIndex(
                infile = this.reference.absoluteFile,
                smrt = true,
                iso = true
        )

//        val results = mutableListOf<String>()


        val bam = SJIndex(
                infile = this.input,
                filter = 0,
                silent = !this.show,
                smrt = true
        )

        GeneReadsCoupler(
                bamIndex = bam,
                reference = ref,
                overlapLevel = this.overlapOfRefReads
        ).saveTo( this.output )

    }
}