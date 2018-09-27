package splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.CliktError
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import splice4k.index.BamIndex
import splice4k.index.GffIndex
import splice4k.index.GtfIndex
import splice4k.index.SJIndex
import splice4k.tools.IdentifyAS
import java.io.FileNotFoundException
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180927
 */



class SGS: CliktCommand(help = "Find AS from NGS") {
    private val input by option(
            "-i",
            "--input",
            help = "Path to input Bam file"
    ).file(exists = true)


    private val reference by option(
            "-r",
            "--reference",
            help = "Path to reference file"
    ).file(exists = true)


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


    private val spliceJunction by option(
            "-j",
            "-junctions",
            help = "Path to extracted Splice junctions file by Splice4k or STAR"
    ).file(exists = true)


    private val isStar by option (
            "-s",
            "--star",
            help = "Is STAR SJ.out.tab file [Only work with -j|--junctions parameter]"
    ).flag(default = false)


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


    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file()


    private val show by option(
            "-v",
            "-verbose",
            help = "Enable detailed messages"
    ).flag(default = true)


    override fun run() {
        if (
                this.input == null && this.spliceJunction == null ||
                this.output == null ||
                this.reference == null
        ) {
            throw FileNotFoundException("Please set input file or output directory")
        }

        if ( this.output == null ) {
            throw CliktError("-o not added")
        }

        val sj = when {
            this.spliceJunction != null -> {
                SJIndex(
                        infile = this.spliceJunction.toString(),
                        filter = this.junctionsFilter,
                        star = this.isStar
                )
            }
            else -> {
                BamIndex(
                        infile = this.input.toString(),
                        filter = this.junctionsFilter
                )
            }
        }

        val ref = when {
            Regex(".*\\.gtf$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                GtfIndex(this.reference.toString())
            }

            Regex(".*\\.gff3?$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                GffIndex(this.reference.toString())
            }
            else -> {
                exitProcess(0)
            }
        }

        val identifyAS = IdentifyAS( overlapOfExonIntron = this.overlapOfExonIntron )

        identifyAS.writeTo(
                outfile = output!!,
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