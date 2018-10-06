package com.splice4k.cli


import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.base.Exons
import com.splice4k.base.SpliceEvent
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.tools.IdentifyAS
import java.io.PrintWriter
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20181006
 */



class SGS: CliktCommand(help = "Find AS from NGS") {
    private val input by argument(
            help = "Path to input file, multiple files separate by space [bam|sam|sj|star SJ.out.tab]"
    ).file(exists = true).multiple()

    private val bam by option(
            "-b",
            "--bam",
            help = "Bam file for calculation of IR PSI value, index needed"
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


    private val unique by option(
            "-u",
            "--unique",
            help = "Unique of STAR"
    ).flag(default = false)


    private val junctionsFilter by option(
            "-c",
            "--count",
            help = "Filter low abundance junctions [default: 3]"
    ).int().default(3).validate {
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

        if ( this.output.isDirectory ) {
            println("Please set path of output file [event not exists]")
            exitProcess(0)
        }

        if ( this.input.isEmpty() ) {
            println("input files are required")
        } else {
            val psis = mutableMapOf<SpliceEvent, MutableMap<String, Double?>>()
            val labels = mutableListOf<String>()
            val results = mutableMapOf<SpliceEvent, MutableList<Exons>>()

            val ref = AnnotationIndex(
                    infile = this.reference.absoluteFile,
                    smrt = false
            )


            for ( it in this.input ) {
                labels.add( it.name )

                val sj = SJIndex(
                        infile = it.absoluteFile,
                        filter = this.junctionsFilter,
                        silent = !this.show,
                        smrt = false,
                        unique = this.unique
                )


                val bamFile = when ( sj.fileFormat == "bam" ) {
                    true -> it
                    else -> this.bam
                }

                val identifyAS = IdentifyAS(
                        overlapOfExonIntron = this.overlapOfExonIntron,
                        bamFile = bamFile
                )

                val data = identifyAS.matchEventsWithRef(
                        event = sj.data.values.toList(),
                        annotations = ref.data,
                        error = this.error,
                        threads = this.threads,
                        show = this.show
                )

                for ( (k, values) in data ) {
                    if ( results.containsKey(k) ) {
                        results[k]!!.addAll(values)
                    } else {
                        results[k] = values.toMutableList()
                    }

                    if ( psis.containsKey(k) ) {
                        psis[k]!![labels.last()] = k.psi
                    } else {
                        psis[k] = mutableMapOf(labels.last() to k.psi)
                    }
                }

            }


            if ( !this.output.absoluteFile.parentFile.exists() ) {
                this.output.absoluteFile.parentFile.mkdirs()
            }

            val writer = PrintWriter(this.output)
            val tmpResults = mutableSetOf<String>()

            writer.println("#spliceRange\tspliceType\tspliceSites\tgene\ttranscript\texon\t${labels.joinToString("\t")}")

            for ((key, values) in results ) {
                for ( v in values ) {
                    val psi = mutableListOf<Double?>()
                    for ( label in labels ) {
                        psi.add(psis[key]!![label])
                    }
                    tmpResults.add("$key\t$v\t${psi.joinToString("\t")}")
                }
            }
            writer.print(tmpResults.asSequence().sorted().distinct().joinToString("\n"))
            writer.close()
        }

    }
}