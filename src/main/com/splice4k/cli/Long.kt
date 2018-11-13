package com.splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.base.SpliceEvent
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.sms.tools.SJFinder
import com.splice4k.sms.tools.TranscriptsReadsCoupler
import com.splice4k.tools.CountTable
import com.splice4k.tools.PSITable
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20181006
 */


class Long: CliktCommand(help = "Identify alternative splicing events from SMRT-seq") {

    private val input by argument(
            help = "Path to input file, multiple files separate by space [BAM|SAM|gmap alignments(-A without -f parameters)|SJ]"
    ).file(exists = true).multiple()


    private val reference by option(
            "-r",
            "--reference",
            help = "Path to reference file [gtf|gff3]"
    ).file(exists = true).required()

    private val bamFile by option(
            "-b",
            "--bam",
            help = "Path to BAM/SAM file"
    ).file(exists = true)


    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file().required()


    private val outputPSITable by option(
            "--psi-table",
            help = "Output PSI table of same starts and same ends junctions"
    ).flag(default = false)

    private val outputCountTable by option(
            "--count-table",
            help = "Output Count table of junctions"
    ).flag(default = false)


    private val threads by option(
            "-p",
            "--process",
            help = "Number of processes to use"
    ).int().default(1).validate {
        require( it > 0 && it <= Runtime.getRuntime().availableProcessors())
    }


    private val error by option(
            "-e",
            help = "The error to identify whether AS event exists [default: 3bp]"
    ).int().default(3).validate {
        require( it >= 0 ) {"this value must be positive"}
    }


    private val junctionsFilter by option(
            "-c",
            "--count",
            help = "Filter low abundance junctions [default: 1]"
    ).int().validate {
        it >= 1
    }

    private val overallJunctionFilter by option(
            "--overall-count",
            help="""
                Filter low abundance junctions across all samples.
                Eg: set this parameter to 100, then junctions that total counts in all samples are lower than 100 will filtered.
                Note: this parameter won't disable -c|--count, when this parameter is used, the -c will set to 1 by default
            """
    ).int().validate {
        it >= 0
    }


    private val overlapOfRefReads by option(
            "--overlap-ref-reads",
            help = "Minimal overlap level to match reference with reads [default: 90.0]"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }


    private val overlapOfExonIntron by option(
            "--overlap-exon-intron",
            help = "Minimal overlap level between exon with intron required for intron retention identification [default: 90.0]"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }


    private val show by option(
            "-v",
            "--verbose",
            help = "Enable detailed messages"
    ).flag(default = false)


    override fun run() {

        if ( this.output.isDirectory ) {
            println("Please set path of output file [event not exists]")
            exitProcess(0)
        }

        if ( this.input.isEmpty() ) {
            println("Input files are required")
            exitProcess(0)
        }

        val logger = Logger.getLogger(Long::class.java)

        // 生成各种文件路径
        if (!this.output.absoluteFile.parentFile.exists()) this.output.absoluteFile.parentFile.mkdirs()

        val psis = mutableMapOf<SpliceEvent, MutableMap<String, String>>()
        val labels = mutableListOf<String>()
        val results = mutableMapOf<SpliceEvent, MutableList<String>>()
        val psiTable = PSITable()
        val subtypes = mutableMapOf<SpliceEvent, HashSet<String>>()
        val junctions = mutableListOf<Pair<SJIndex, File?>>()

        val ref = AnnotationIndex(
                infile = this.reference.absoluteFile,
                smrt = true
        )

        for ( it in this.input ) {
            labels.add(it.name)

            val sj =  SJIndex(
                    infile = it.absoluteFile,
                    silent = !this.show,
                    smrt = true,
                    filter = this.junctionsFilter ?: 0
            )

            if ( this.outputPSITable ) {
                psiTable.addJunctionGraph(sj)
            }

            junctions.add(
                    Pair(
                            sj,
                            when( sj.fileFormat ) {
                                "bam" -> sj.infile
                                else -> this.bamFile
                            }
                    )
            )
        }

        if ( this.outputCountTable ) {
            val countTable = CountTable()

            countTable.writeTo(
                    this.output,
                    junctions.map { it.first }
            )
        }

        var overallFiltered: HashSet<String>? = null
        if ( this.overallJunctionFilter != null ) {
            val countTable = CountTable()
            overallFiltered = countTable.filter(
                    junctions.map { it.first },
                    this.overallJunctionFilter!!
            )
        }

        logger.info("Start to compare ref and reads")
        for ( (sj, bamFile) in junctions ) {

            val matched = TranscriptsReadsCoupler(
                    reference = ref,
                    reads = sj,
                    overlap = this.overlapOfRefReads,
                    distanceError = this.error
            )

            logger.info("Start to identify splice events")
            logger.info("Predicting Alternative Splicing events of ${sj.infile.name}")
            val data = SJFinder(
                    template = matched,
                    bamIndex = sj,
                    refIndex = ref,
                    silent = !this.show,
                    overlapOfExonIntron = this.overlapOfExonIntron,
                    error = this.error,
                    threads = this.threads,
                    bamFile = bamFile,
                    overallFiltered = overallFiltered
            ).results

            labels.add( sj.infile.name )
            for ((k, values) in data) {
                val tmpRes = results[k] ?: mutableListOf()
                tmpRes.addAll(values)
                results[k] = tmpRes

                val tmpPSI = psis[k] ?: mutableMapOf()
                tmpPSI[labels.last()] = k.getPsi()
                psis[k] = tmpPSI

                val tmpSubtypes = subtypes[k] ?: hashSetOf()
                if ( k.subtypes != "NA" ) {
                    tmpSubtypes.add(k.subtypes)
                }
                subtypes[k] = tmpSubtypes
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
            for ( v in values ) {
                val psi = mutableListOf<String>()
                for ( label in labels ) {
                    psi.add(psis[key]!![label] ?: "NA")
                }

                var tmpSubtype = subtypes[key]?.joinToString(separator = ",") ?: "NA"
                if ( tmpSubtype == "" ) {
                    tmpSubtype = "NA"
                }

                tmpResults.add(
                        "${key.getString(tmpSubtype)}\t" +
                        "${if (key.isNovel) 1 else 0}\t" +
                        "$v\t" +
                        psi.joinToString("\t")
                )
            }
        }
        writer.print(tmpResults.asSequence().sorted().distinct().joinToString("\n"))
        writer.close()

        if ( this.outputPSITable ) {
            psiTable.writeTo( this.output.absoluteFile )
        }
    }
}
