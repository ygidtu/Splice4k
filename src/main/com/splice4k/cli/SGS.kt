package com.splice4k.cli


import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.tools.IdentifyAS


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180929
 */



class SGS: CliktCommand(help = "Find AS from NGS") {
    private val input by option(
            "-i",
            "--input",
            help = "Path to input file [bam|sam|sj|star SJ.out.tab]"
    ).file(exists = true).required()

    private val bam by option(
            "-b",
            "--bam",
            help = "Bam file for calculation of IR PSI value"
    ).file(exists = true)


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

        val sj = SJIndex(
                infile = this.input.absoluteFile,
                filter = this.junctionsFilter,
                silent = !this.show,
                smrt = false
        )


        val ref = AnnotationIndex(
                infile = this.reference.absoluteFile,
                smrt = false
        )

        val bamFile = when ( sj.fileFormat == "bam" ) {
            true -> this.input
            else -> this.bam
        }

        val identifyAS = IdentifyAS(
                overlapOfExonIntron = this.overlapOfExonIntron,
                bamFile = bamFile
        )

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