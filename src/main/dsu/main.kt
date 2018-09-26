package dsu


import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.CliktError
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import dsu.second.identifier.IdentifyAS
import dsu.second.index.BamIndex
import dsu.second.index.GffIndex
import dsu.second.index.GtfIndex
import dsu.second.index.SJIndex
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
            println("20180925")
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
    private val output by option("-o", "--output", help = "Path to output file").file()

    private val error by option(
            "-e",
            help = "Chromosome coordinate error"
    ).int().default(3).validate {
        require( it >= 0 ) {"this value must be positive"}
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
        if (!this.output!!.absoluteFile.parentFile.exists()) this.output!!.absoluteFile.parentFile.mkdirs()

        val bamFile = this.input!!.absoluteFile

        logger.info("Start to read $bamFile")
        val bam = BamExtractor(bamFile.toString(), silent = !this.show)

        val ref : Extractor

        logger.info("Start to read ${this.reference}")
        val refFile = this.reference!!.absoluteFile
        // val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")
        when {

            Regex(".*\\.gff3?$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GffExtractor(refFile.toString(), !this.show)
            }

            Regex(".*\\.gtf$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GtfExtractor(refFile.toString(), !this.show)
            }

            Regex(".*\\.(bam|sam)$").matches(
                    this.reference.toString().toLowerCase()
            ) -> {
                ref = GtfExtractor(refFile.toString(), !this.show)
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
                error = this.error,
                overlapOfExonIntron = this.overlapOfExonIntron
        ).saveTo(this.output!!.absoluteFile.toString())
    }
}



class Short: CliktCommand(help = "Find AS from NGS") {
    val input by option("-i", "--input", help = "Path to input Bam file").file(exists = true)

    val reference by option("-r", "--reference", help = "Path to reference file")

    val spliceJunction by option(
            "-j",
            "-junctions",
            help = "Path to extracted Spolice junctions file"
    ).file(exists = true)



    private val overlapOfExonIntron by option(
            "--overlap-exon-intron",
            help = "Minimal overlap level between exon with intron required for intron retention identification [default: 90.0]"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }

    private val error by option(
            "-e",
            help = "Chromosome coordinate error [default: 3]"
    ).int().default(3).validate {
        require( it >= 0 ) {"this value must be positive"}
    }


    private val output by option("-o", "--output", help = "Path to output file").file()
    private val show by option("--show", "-s", help = "Enable detailed messages").flag(default = false)

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
                SJIndex(this.spliceJunction.toString())
            }
            else -> {
                BamIndex(this.input.toString())
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

        val identifyAS = IdentifyAS(
                overlapOfExonIntron = this.overlapOfExonIntron,
                silent = !this.show,
                distanceError = this.error
        )

        identifyAS.writeTo(
                outfile = output!!,
                results = identifyAS.matchEventsWithRef(event = sj, annotation = ref)
        )

    }
}


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters().subcommands(Extract()).subcommands(dsu.Short()).subcommands(dsu.Long())
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
