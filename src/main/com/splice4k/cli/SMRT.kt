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
import com.splice4k.smrt.tools.GeneReadsCoupler
import com.splice4k.smrt.tools.SJFinder
import org.apache.log4j.FileAppender
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.apache.log4j.PatternLayout
import java.io.File
import java.io.PrintWriter
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180928
 */



class SMRT: CliktCommand(help = "Find AS from PacBio data") {

    private val input by argument(
            help = "Path to input Bam/Sam file"
    ).file(exists = true).multiple()


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
            help = "Filter low abundance junctions [default: 0]"
    ).int().default(0).validate {
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


    private val log by option(help = "Path to the log file").file()


    /**
     * 为log添加文件appender
     * @param logFile log文件的地址
     */
    private fun addFileAppender(logFile: File) {
        val fa = FileAppender()
        fa.name = "FileLogger"
        fa.file = logFile.absolutePath
        fa.layout = PatternLayout("[%d{yyyy-MM-dd HH:mm:ss}] [%-5p] [%c{1}:%L] - %m%n")
        fa.threshold = Level.DEBUG
        fa.append = true
        fa.activateOptions()

        Logger.getRootLogger().addAppender(fa)
    }

    override fun run() {
        this.log?.let {
            this.addFileAppender(it)
        }

        if ( this.output.isDirectory ) {
            println("Please set path of output file [event not exists]")
            exitProcess(0)
        }

        val logger = Logger.getLogger(SMRT::class.java)

        // 生成各种文件路径
        if (!this.output.absoluteFile.parentFile.exists()) this.output.absoluteFile.parentFile.mkdirs()

        val psis = mutableMapOf<SpliceEvent, MutableMap<String, Double?>>()
        val labels = mutableListOf<String>()
        val results = mutableMapOf<SpliceEvent, MutableList<String>>()


        // val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")
        val ref = AnnotationIndex(
                infile = this.reference.absoluteFile,
                smrt = true
        )


        for ( it in this.input ) {
            val bam =  SJIndex(
                    infile = it.absoluteFile,
                    silent = !this.show,
                    smrt = true,
                    filter = this.junctionsFilter
            )

            logger.info("Start to compare ref and reads")
            val matched = GeneReadsCoupler(
                    reference = ref,
                    reads = bam,
                    overlap = this.overlapOfRefReads,
                    distanceError = this.error
            )

            logger.info("Start to identify splice events")
            val data = SJFinder(
                    template = matched,
                    bamIndex = bam,
                    refIndex = ref,
                    silent = !this.show,
                    overlapOfExonIntron = this.overlapOfExonIntron,
                    error = this.error,
                    threads = this.threads
            ).results

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

            if ( !this.output.parentFile.exists() ) {
                this.output.parentFile.mkdirs()
            }

            val writer = PrintWriter(this.output)

            writer.println("#spliceRange\tspliceType\tspliceSites\tgene\ttranscript\texon\t${labels.joinToString("\t")}")

            for ((key, values) in results ) {
                for ( v in values ) {
                    val psi = mutableListOf<Double?>()
                    for ( label in labels ) {
                        psi.add(psis[key]!![label])
                    }
                    writer.println("$key\t$v\t${psi.joinToString("\t")}")
                }
            }
        }
    }
}
