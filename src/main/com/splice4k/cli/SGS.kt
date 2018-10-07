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
import me.tongfei.progressbar.ProgressBar
import java.io.PrintWriter
import java.io.File
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
            help = "Path to BAM/SAM files, or the directory contains BAM/SAM files. [default: current running directory]\n" +
                    "\tIf input files are BAM/SAM files, this parameter won't work\n" +
                    "\tIf specified path to directory contains BAM/SAM files corresponding to STAR SJ.out.tab files, this program will auto match those files\n" +
                    "\tIf specified BAM/SAM file with this parameter, then this program will calculate PSI of IR using this file\n"
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
                        smrt = false
                )


                val bamFile = when ( sj.fileFormat ) {
                    "bam" -> it
                    "star" -> {

                        if ( this.bam != null ) {
                            if ( this.bam!!.isFile ) {   // bam is file, this use this file
                                this.bam
                            } else {                     // this.bam is directory, then try to find the corresponding bam file
                                var bamFile = this.bam

                                val pattern = ".*${it.name.replace("[_\\.]SJ.out.tab".toRegex(), "")}[\\._](\\w+\\.)*bam$"
                                        .toRegex(RegexOption.IGNORE_CASE)

                                for ( i in this.bam!!.walkTopDown() ) {
                                    if ( i.name.matches( pattern ) ) {
                                        bamFile = i.absoluteFile
                                        break
                                    }
                                }
                                bamFile
                            }
                        } else {
                            this.bam
                        }
                    }
                    else -> this.bam
                }

                if ( sj.fileFormat == "star" && bamFile != null ) {
                    println( "${sj.infile.name} -> ${bamFile.name}" )
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

                for ( (k, values) in ProgressBar.wrap(data.iterator(), "") ) {
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