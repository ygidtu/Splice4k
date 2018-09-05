package dsu

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int


import dsu.second.index.BamIndex
import dsu.third.extractor.BamExtractor
import dsu.third.extractor.Extractor
import dsu.third.extractor.GffExtractor
import dsu.third.extractor.GtfExtractor
import dsu.third.template.GeneReadsCoupler
import dsu.third.template.SJFinder
import org.apache.log4j.FileAppender
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.apache.log4j.PatternLayout
import java.io.File
import kotlin.system.exitProcess



class Parameters: CliktCommand() {
    private val version by option("--version", "-v", help = "version").flag(default = false)

    override fun run() {
        if ( this.version ) {
            exitProcess(0)
        }
    }
}


class Extract: CliktCommand(help = "Extract junctions from Bam/Sam files") {
    private val input by argument("-i", help = "Path to input Bam/Sam file").file(exists = true)
    private val output by argument("-o", help = "Path to output file").file()

    override fun run() {
        BamIndex(this.input.toString()).writeTo(this.output)
    }
}


class PacBio: CliktCommand(help = "Find AS from PacBio data") {

    private val input by argument("-i", help = "Path to input Bam/Sam file").file(exists = true)
    private val reference by argument("-r", help = "Path to reference file [gtf|gff3]").file(exists = true)
    private val output by argument("-o", help = "Path to output directory").file()

    private val error by option(
            "-e",
            help = "Chromosome coordinate error"
    ).int().default(3).validate {
        require( it >= 0 ) {"this value must be positive"}
    }

    private val foldChange by option(
            "--fold-change",
            "-fc",
            help = "Minimal fold change required to match genes with same reads [default: 1.5]"
    ).double().default(1.5).validate {
        require( it > 0 ) {"this value must be positive"}
    }

    private val overlapOfRefReads by option(
            "--overlap-ref-reads",
            help = "Minimal overlap level to match reference with reads"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }

    private val overlapOfExonIntron by option(
            "--overlap-exon-intron",
            help = "Minimal overlap level between exon with intron required for intron retention identification"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }

    private val silent by option("--show", "-s", help = "Enable detailed messages").flag(default = false)
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


    override fun run() {
        if ( this.log != null ) {
            this.addFileAppender(log.toString())
        }

        val logger = Logger.getLogger(PacBio::class.java)

        // 生成各种文件路径
        val outDir = this.output.absoluteFile

        if (!outDir.exists()) outDir.parentFile.mkdirs()

        val bamFile = this.input.absoluteFile
        val bamTsv = File(outDir, bamFile.name.split(".")[0] + ".tsv")

        logger.info("Start to read $bamFile")
        val bam = BamExtractor(bamFile.toString(), silent = silent)

        bam.saveTo(bamTsv.toString())

        val ref : Extractor

        logger.info("Start to read ${this.reference}")
        val refFile = this.reference.absoluteFile
        val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")
        when {

            Regex(".*\\.gff3?$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GffExtractor(refFile.toString(), silent)
                ref.saveTo(refTsv.toString())
            }

            Regex(".*\\.gtf$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GtfExtractor(refFile.toString(), silent)
                ref.saveTo(refTsv.toString())
            }

            Regex(".*\\.(bam|sam)$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GtfExtractor(refFile.toString(), silent)
                ref.saveTo(refTsv.toString())
            }

            else -> {
                println("Input should be gtf or gff3 format")
                exitProcess(2)
            }

        }

        logger.info("Start to compare ref and reads")
        val matched = GeneReadsCoupler(
                ref, bam,
                overlap = this.overlapOfRefReads,
                foldChange = this.foldChange,
                distanceError = this.error
        )

        matched.saveTo(File(outDir, "gene_reads_pairs.tsv").toString())
        matched.saveTemplate(File(outDir, "templates.tsv").toString())
        matched.saveNovel(File(outDir, "novel.tsv").toString())
        matched.savePotentialFusion(File(outDir, "potential_fusions.tsv").toString())

        logger.info("Start to identify splice events")
        SJFinder(
                matched,
                distance = this.error,
                overlap = this.overlapOfExonIntron
        ).saveTo(File(outDir, "final.tsv").toString())
    }
}


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters().subcommands(Extract()).subcommands(PacBio())
    if (args.size <= 1) {
        val help = when {
            args.isEmpty() -> cmd.getFormattedHelp()
            args[0] == "pacbio" -> PacBio().getFormattedHelp()
            args[0] == "extract" -> Extract().getFormattedHelp()
            else -> cmd.getFormattedHelp()
        }
        println(help)
        exitProcess(0)
    }

    try {
        cmd.parse(args)
    } catch (e: CliktError) {
        logger.error(e.localizedMessage)
        println(cmd.getFormattedHelp())
    }
}