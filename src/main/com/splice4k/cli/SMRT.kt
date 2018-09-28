package com.splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.index.BamIndex
import com.splice4k.index.GffIndex
import com.splice4k.index.GtfIndex
import com.splice4k.smrt.tools.GeneReadsCoupler
import com.splice4k.smrt.tools.SJFinder
import org.apache.log4j.FileAppender
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.apache.log4j.PatternLayout
import java.io.File
import java.io.FileNotFoundException
import kotlin.system.exitProcess
import com.splice4k.tools.FileValidator


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180928
 */



class SMRT: CliktCommand(help = "Find AS from PacBio data") {

    private val input by option(
            "-i",
            "--input",
            help = "Path to input Bam/Sam file"
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

        val logger = Logger.getLogger(SMRT::class.java)
        val fileValidator = FileValidator()

        // 生成各种文件路径
        if (!this.output.absoluteFile.parentFile.exists()) this.output.absoluteFile.parentFile.mkdirs()

        val bam = when(fileValidator.check(this.input)) {
            "bam" -> BamIndex(
                    infile = this.input.absoluteFile.toString(),
                    silent = !this.show,
                    smrt = true,
                    filter = this.junctionsFilter
            )
            else -> {
                logger.info("Please check input file format")
                exitProcess(2)
            }
        }

        // val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")
        val ref = when(fileValidator.check(this.reference)) {

             "gff" -> GffIndex(
                        infile = this.reference.absoluteFile.toString(),
                        smrt = true
                )

             "gtf" -> GtfIndex(
                        infile = this.reference.absoluteFile.toString(),
                        smrt = true
                )

            else -> {
                logger.info("Please check reference file format")
                exitProcess(2)
            }

        }

        logger.info("Start to compare ref and reads")
        val matched = GeneReadsCoupler(
                reference = ref,
                reads = bam,
                overlap = this.overlapOfRefReads,
                distanceError = this.error
        )

        logger.info("Start to identify splice events")
        SJFinder(
                template = matched,
                bamIndex = bam,
                refIndex = ref,
                silent = !this.show,
                overlapOfExonIntron = this.overlapOfExonIntron,
                error = this.error,
                threads = this.threads
        ).saveTo(this.output.absoluteFile.toString())
    }
}
