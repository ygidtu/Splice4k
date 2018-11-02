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
import com.splice4k.tools.PSITable
import org.apache.log4j.Logger
import java.io.PrintWriter
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20181006
 */


class SGS: CliktCommand(help = "Identify alternative splicing events from RNA-seq") {
    private val input by argument(
            help = "Path to input file, multiple files separate by space [BAM|SAM|STAR SJ.out.tab|gmap align|SJ]"
    ).file(exists = true).multiple()

    private val bam by option(
            "-b",
            "--bam",
            help = """
                Path to BAM/SAM files, or the directory contains BAM/SAM files. [default: current running directory]
                - If input files are BAM/SAM files, this parameter won't work
                - If specified path to directory contains BAM/SAM files corresponding to STAR SJ.out.tab files, this program will auto match those files
                - If specified BAM/SAM file with this parameter, then this program will calculate PSI of IR using this file
                """
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

    private val outputPSITable by option(
            "--psi",
            help = "Output PSI table of same starts and same ends junctions"
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
        val logger = Logger.getLogger(SGS::class.java)
        if ( this.output.isDirectory ) {
            logger.error("Please set path of output file [event not exists]")
            exitProcess(0)
        }

        if ( this.input.isEmpty() ) {
            logger.error("input files are required")
        } else {
            val psis = mutableMapOf<SpliceEvent, MutableMap<String, String>>()
            val labels = mutableListOf<String>()
            val results = mutableMapOf<SpliceEvent, MutableList<Exons>>()
            val psiTable = PSITable()

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


                if ( this.outputPSITable ) {
                    psiTable.addJunctionGraph(sj)
                }


                val bamFile = when ( sj.fileFormat ) {
                    "bam" -> it
                    "star" -> {

                        val bamDirectory = this.bam ?: it.absoluteFile.parentFile

                        if ( bamDirectory!!.isFile ) {   // bam is file, this use this file
                            bamDirectory
                        } else {                     // this.bam is directory, then try to find the corresponding bam file
                            var bamFile = bamDirectory

                            val pattern = ".*${it.name.replace("[_.]?SJ.out.tab".toRegex(), "")}[._]?(\\w+.)*bam$"
                                    .toRegex(RegexOption.IGNORE_CASE)

                            for ( i in bamDirectory.walkTopDown() ) {
                                if ( i.name.matches( pattern ) ) {
                                    bamFile = i.absoluteFile
                                    break
                                }
                            }
                            bamFile
                        }
                    }
                    else -> this.bam
                }

                if ( sj.fileFormat == "star" && bamFile != null ) {
                    logger.info( "${sj.infile.name} -> ${bamFile.name}" )
                }

                val identifyAS = IdentifyAS(
                        overlapOfExonIntron = this.overlapOfExonIntron,
                        bamFile = bamFile
                )

                val data = identifyAS.matchEventsWithRef(
                        event = sj.data.values.toList(),
                        annotations = ref,
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
                        psis[k]!![labels.last()] = k.getPsi()
                    } else {
                        psis[k] = mutableMapOf(labels.last() to k.getPsi())
                    }
                }
            }


            if ( !this.output.absoluteFile.parentFile.exists() ) {
                this.output.absoluteFile.parentFile.mkdirs()
            }

            val writer = PrintWriter(this.output)
            val tmpResults = mutableSetOf<String>()

            writer.println(
                    "#spliceRange\t" +
                    "spliceType\t" +
                    "subtype\t" +
                    "spliceSites\t" +
                    "isNovel\t" +
                    "gene\t" +
                    "transcript\t" +
                    "exon\t" +
                    labels.joinToString("\t")
            )

            for ((key, values) in results ) {
                val psi = mutableListOf<String>()
                for ( label in labels ) {
                    psi.add(psis[key]!![label] ?: "NA")
                }

                val gene = mutableSetOf<String>()
                val transcript = mutableSetOf<String>()
                val exon = mutableSetOf<String>()

                for (v in values) {
                    v.source["gene"]?.let {
                        gene.addAll(it)
                    }

                    v.source["transcript"]?.let {
                        transcript.addAll(it)
                    }

                    exon.add(v.exonId)
                }

                if ( gene.isEmpty() ) gene.add("NA")
                if ( transcript.isEmpty() ) transcript.add("NA")
                if ( exon.isEmpty() ) exon.add("NA")

                tmpResults.add(
                        "$key\t" +
                        "${if (key.isNovel) 1 else 0}\t" +
                        "${gene.joinToString(",")}\t" +
                        "${transcript.joinToString(",")}\t" +
                        "${exon.joinToString(",")}\t" +
                        psi.joinToString("\t")
                )

            }
            writer.print(tmpResults.asSequence().sorted().distinct().joinToString("\n"))
            writer.close()


            if ( this.outputPSITable ) {
                psiTable.writeTo( this.output.absoluteFile )
            }
        }
    }
}