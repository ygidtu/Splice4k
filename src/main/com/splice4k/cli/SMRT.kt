package splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import org.apache.log4j.FileAppender
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.apache.log4j.PatternLayout
import splice4k.index.BamIndex
import splice4k.index.GffIndex
import splice4k.index.GtfIndex
import splice4k.smrt.tools.GeneReadsCoupler
import splice4k.smrt.tools.SJFinder
import java.io.File
import java.io.FileNotFoundException
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180927
 */



class SMRT: CliktCommand(help = "Find AS from PacBio data") {

    private val input by option(
            "-i",
            "--input",
            help = "Path to input Bam/Sam file"
    ).file(exists = true)


    private val reference by option(
            "-r",
            "--reference",
            help = "Path to reference file [gtf|gff3]"
    ).file(exists = true)


    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file()


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
            help = "Filter low abundance junctions [default: 10]"
    ).int().default(10).validate {
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


    private val log: String? by option(help = "Path to the log file")


    /**
     * 为log添加文件appender
     * @param logFile log文件的地址
     */
    private fun addFileAppender(logFile: String) {
        val fa = FileAppender()
        fa.name = "FileLogger"
        fa.file = File(logFile).absolutePath
        fa.layout = PatternLayout("[%d{yyyy-MM-dd HH:mm:ss}] [%-5p] [%c{1}:%L] - %m%n")
        fa.threshold = Level.DEBUG
        fa.append = true
        fa.activateOptions()

        Logger.getRootLogger().addAppender(fa)
    }


    private fun checkFile(infile: File?) {
        if ( infile == null ) {
            throw FileNotFoundException("Input could not be null")
        }
    }

    override fun run() {
        if ( this.log != null ) {
            this.addFileAppender(log.toString())
        }

        this.checkFile(this.input)
        this.checkFile(this.output)
        this.checkFile(this.reference)

        val logger = Logger.getLogger(Long::class.java)

        // 生成各种文件路径
        if (!this.output!!.absoluteFile.parentFile.exists()) this.output!!.absoluteFile.parentFile.mkdirs()

        val bamFile = this.input!!.absoluteFile

        logger.info("Start to read $bamFile")
        val bam = BamIndex(
                infile = bamFile.toString(),
                silent = !this.show,
                smrt = true,
                filter = this.junctionsFilter
        )

        logger.info("Start to read ${this.reference}")
        val refFile = this.reference!!.absoluteFile
        // val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")
        val ref = when {

            Regex(".*\\.gff3?$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                GffIndex(
                        infile = refFile.toString(),
                        smrt = true
                )
            }

            Regex(".*\\.gtf$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                GtfIndex(
                        infile = refFile.toString(),
                        smrt = true
                )
            }

            else -> {
                println("Input should be gtf or gff3 format")
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
        ).saveTo(this.output!!.absoluteFile.toString())
    }
}
