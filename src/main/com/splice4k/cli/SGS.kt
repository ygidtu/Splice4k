package com.splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
//import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.index.BamIndex
import com.splice4k.index.GffIndex
import com.splice4k.index.GtfIndex
import com.splice4k.index.SJIndex
import com.splice4k.tools.IdentifyAS
import kotlin.system.exitProcess
import com.splice4k.tools.FileValidator
import org.apache.log4j.Logger


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180928
 */



class SGS: CliktCommand(help = "Find AS from NGS") {
    private val input by option(
            "-i",
            "--input",
            help = "Path to input Bam file"
    ).file(exists = true).required()


    private val reference by option(
            "-r",
            "--reference",
            help = "Path to reference file"
    ).file(exists = true).required()


    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file().required()


//    private val format by option(
//            "-f",
//            "--format",
//            help = "Input file format [bam|sj|star].\nbam -> BAM/SAM file \nsj -> extracted splice junctions by this program.\nstar -> SJ.out.tab file from STAR"
//    ).choice("bam", "sj", "star").default("bam")


    private val junctionsFilter by option(
            "-c",
            "--count",
            help = "Filter low abundance junctions [default: 5]"
    ).int().default(5).validate {
        it >= 0
    }


    private val error by option(
            "-e",
            help = "The error to identify whether AS event exists [default: 3bp]"
    ).int().default(3).validate {
        require( it >= 0 ) {"this value must be positive"}
    }


    private val threads by option(
            "-p",
            "--process",
            help = "Number of processes to use"
    ).int().default(1).validate {
        require( it > 0 && it <= Runtime.getRuntime().availableProcessors())
    }


    private val overlapOfExonIntron by option(
            "--overlap-exon-intron",
            help = "Minimal overlap level between exon with intron required for intron retention identification [default: 90.0]"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }


    private val show by option(
            "-v",
            "-verbose",
            help = "Enable detailed messages"
    ).flag(default = false)


    override fun run() {

        val logger = Logger.getLogger(SGS::class.java)
        val fileValidator = FileValidator()
        var format = fileValidator.check(this.input)

        val sj = when(format) {
            "sj" -> SJIndex(
                    infile = this.input.absolutePath.toString(),
                    filter = this.junctionsFilter,
                    star = false
            )

            "bam" -> BamIndex(
                    infile = this.input.absolutePath.toString(),
                    filter = this.junctionsFilter
            )

            "star" -> SJIndex(
                    infile = this.input.absolutePath.toString(),
                    filter = this.junctionsFilter,
                    star = true
            )

            else -> {
                logger.info("Please check the input file format")
                exitProcess(1)
            }
        }

        format = fileValidator.check(this.reference)
        val ref = when(format) {
            "gtf" ->  GtfIndex(this.reference.absolutePath.toString())

            "gff" -> GffIndex(this.reference.absolutePath.toString())

            else -> {
                logger.info("Please check reference file format")
                exitProcess(0)
            }
        }

        val identifyAS = IdentifyAS( overlapOfExonIntron = this.overlapOfExonIntron )

        identifyAS.writeTo(
                outfile = this.output,
                results = identifyAS.matchEventsWithRef(
                        event = sj.data.values.toList(),
                        annotations = ref.data,
                        error = this.error,
                        threads = this.threads,
                        show = this.show
                )
        )
    }
}