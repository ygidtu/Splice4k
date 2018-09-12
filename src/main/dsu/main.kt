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


import dsu.second.identifier.IdentifyAS
import dsu.second.index.*
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
import java.io.FileNotFoundException
import kotlin.system.exitProcess



class Parameters: CliktCommand(invokeWithoutSubcommand = true) {
    private val version by option("--version", "-v", help = "version").flag(default = false)

    override fun run() {
        if ( this.version ) {
            println("2080905")
            exitProcess(0)
        }
    }
}


class Extract: CliktCommand(help = "Extract junctions from Bam/Sam files") {
    private val input by argument(help = "Path to input Bam/Sam file").file(exists = true)
    private val output by argument(help = "Path to output file").file()

    override fun run() {
        BamIndex(this.input.toString()).writeTo(this.output)
    }
}


class Long: CliktCommand(help = "Find AS from PacBio data") {

    private val input by option("-i", "--input", help = "Path to input Bam/Sam file").file(exists = true)
    private val reference by option("-r", "--reference", help = "Path to reference file [gtf|gff3]").file(exists = true)
    private val output by option("-o", "--output", help = "Path to output directory").file()

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

    private val show by option("--show", "-s", help = "Enable detailed messages").flag(default = false)
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
        val outDir = this.output!!.absoluteFile

        if (!outDir.exists()) outDir.parentFile.mkdirs()

        val bamFile = this.input!!.absoluteFile
        val bamTsv = File(outDir, bamFile.name.split(".")[0] + ".tsv")

        logger.info("Start to read $bamFile")
        val bam = BamExtractor(bamFile.toString(), silent = !this.show)

        bam.saveTo(bamTsv.toString())

        val ref : Extractor

        logger.info("Start to read ${this.reference}")
        val refFile = this.reference!!.absoluteFile
        val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")
        when {

            Regex(".*\\.gff3?$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GffExtractor(refFile.toString(), !this.show)
                ref.saveTo(refTsv.toString())
            }

            Regex(".*\\.gtf$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GtfExtractor(refFile.toString(), !this.show)
                ref.saveTo(refTsv.toString())
            }

            Regex(".*\\.(bam|sam)$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GtfExtractor(refFile.toString(), !this.show)
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


class Short: CliktCommand(help = "Find AS from NGS") {
    val input by option("-i", "--input", help = "Path to input Bam file").file(exists = true)
    val reference by option("-r", "--reference", help = "Path to reference file")
    val spliceJunction by option("-j", "-junctions", help = "Path to extracted Spolice junctions file").file(exists = true)
    val output by option("-o", "--output", help = "Path to output file").file()

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

        val sj: SJIndex
        when {
            this.spliceJunction != null -> {
                sj = SJIndex(this.spliceJunction.toString())
            }
            else -> {
                sj = BamIndex(this.input.toString())
                sj.writeTo(File(this.output, "splice_junctions.txt"))
            }
        }

        val ref = when {
            this.reference!!.endsWith(".gtf") -> {
                GtfIndex(this.reference!!)
            }
            this.reference!!.endsWith(".gff") -> {
                GffIndex(this.reference!!)
            }
            else -> {
                exitProcess(0)
            }
        }

        val events = IdentifyAS(
                reference = ref,
                junctions = sj
        )

        events.writeTo(this.output!!)

    }
}


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters().subcommands(Extract()).subcommands(dsu.Long()).subcommands(dsu.Short())
    if (args.size <= 1) {
        val help = when {
            args.isEmpty() -> cmd.getFormattedHelp()
            args[0] == "long" -> dsu.Long().getFormattedHelp()
            args[0] == "extract" -> Extract().getFormattedHelp()
            args[0] == "short" -> dsu.Short().getFormattedHelp()
            else -> cmd.getFormattedHelp()
        }
        println(help)
        exitProcess(0)
    }

    try{
        cmd.parse(args)
        logger.info("Done")
    } catch (e: CliktError) {
        logger.error(e.localizedMessage)
        println(cmd.getFormattedHelp())
        println()
        println(dsu.Long().getFormattedHelp())
        println()
        println(Extract().getFormattedHelp())
    }
}
