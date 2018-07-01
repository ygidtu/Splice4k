package main

import java.io.File
import java.io.IOException

import kotlin.system.exitProcess

import org.apache.commons.cli.Options
import org.apache.commons.cli.DefaultParser
import org.apache.commons.cli.HelpFormatter
import org.apache.commons.cli.ParseException

import org.apache.log4j.FileAppender
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.apache.log4j.PatternLayout

import main.extractor.*
import main.template.*


const val version = "2018.06.28"



/**
 * 为log添加文件appender
 * @param logFile log文件的地址
 */
fun addFileAppender(logFile: String) {
    val fa = FileAppender()
    fa.name = "FileLogger"
    fa.file = File(logFile).absolutePath
    fa.layout = PatternLayout("[%d{yyyy-MM-dd HH:mm:ss}] [%-5p] [%c{1}:%L] - %m%n")
    fa.threshold = Level.DEBUG
    fa.append = true
    fa.activateOptions()

    Logger.getRootLogger().addAppender(fa)
}


/**
 * 设定命令行参数
 */
fun setOptions(): Options {
    val options = Options()

    options.addOption(
            "b",
            "bam",
            true,
            "input [bam|sam] file"
    )
    options.addOption(
            "r",
            true,
            "reference file [gff3|gtf|bam|sam]"
    )
    options.addOption(
            "o",
            "output",
            true,
            "output directory"
    )
    options.addOption(
            "h",
            "help",
            false,
            "print help message"
    )
    options.addOption(
            "s",
            "silent",
            false,
            "suppress output message"
    )
    options.addOption(
            null,
            "log",
            false,
            "write log to file under output directory"
    )
    options.addOption(
            null,
            "fold-change",
            true,
            "Minimal fold change to identify the best match between multiple genes with same reads [default: 1.5]"
    )
    options.addOption(
            null,
            "min-ref-read",
            true,
            "Minimal overlap level to match ref and reads [default: 90.0]"
    )
    options.addOption(
            null,
            "min-AS-length",
            true,
            "Minimal gap between two sites that used to identify alternative splicing [default: 3]"
    )
    options.addOption(
            null,
            "min-exon-intron",
            true,
            "Minimal overlap level to identify the intron retention [default: 90.0]"
    )

    options.addOption(
            "v",
            "version",
            false,
            "display version [current: $version]"
    )

    return options
}


/**
 * @author zhangyiming
 * @since 2018.06.23
 * @version 0.1
 *
 * command line parameters
 */

fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")

    val options = setOptions()

    val parser = DefaultParser()
    val help = HelpFormatter()

    if (args.isEmpty()) {
        help.printHelp("usage message", options)
        exitProcess(0)
    }

    try{
        
        val parm = parser.parse(options, args)

        // add file appender to log
        if (parm.hasOption("log")) {
            addFileAppender(
                    parm.getOptionValue(
                            "log",
                            File(
                                    parm.getOptionValue("output"),
                                    "Splice4k.log"
                            ).toString()
                    )
            )

        }


        var silent = false
        when {
            parm.hasOption("help") -> {
                help.printHelp("usage message", options)
                exitProcess(0)
            }

            parm.hasOption("version") -> {
                println("Current version: $version")
                exitProcess(0)
            }

            parm.hasOption("silent") -> silent = true
        }


        val outDir = File(parm.getOptionValue("output")).absoluteFile

        if (!outDir.exists()) outDir.parentFile.mkdirs()

        val bamFile = File(parm.getOptionValue("bam")).absoluteFile
        val bamTsv = File(outDir, bamFile.name.split(".")[0] + ".tsv")

        println()
        logger.info("Start to read $bamFile")
        val bam = BamExtractor(bamFile.toString(), silent = silent)
        bam.saveTo(bamTsv.toString())


        val ref : Extractor

        println()
        logger.info("Start to read ${parm.getOptionValue("r")}")
        when {

             Regex(".*\\.gff3?$").matches(
                     parm.getOptionValue("r").toLowerCase()
             ) -> {
                val refFile = File(parm.getOptionValue("r")).absoluteFile
                val refTsv = File(outDir, refFile.name.split(".")[0] + ".tsv")

                ref = GffExtractor(refFile.toString(), silent)
                ref.saveTo(refTsv.toString())
            }

            Regex(".*\\.gtf$").matches(
                    parm.getOptionValue("r").toLowerCase()
            ) -> {
                val refFile = File(parm.getOptionValue("r")).absoluteFile
                val refTsv = File(outDir, refFile.name.split(".")[0] + "tsv")

                ref = GtfExtractor(refFile.toString(), silent)
                ref.saveTo(refTsv.toString())
            }

            Regex(".*\\.(bam|sam)$").matches(
                    parm.getOptionValue("r").toLowerCase()
            ) -> {
                val refFile = File(parm.getOptionValue("r")).absoluteFile
                val refTsv = File(outDir, refFile.name.split(".")[0] + "tsv")

                ref = GtfExtractor(refFile.toString(), silent)
                ref.saveTo(refTsv.toString())
            }

            else -> {
                println("Input should be gtf or gff3 format")
                exitProcess(2)
            }

        }

        println()
        logger.info("Start to compare ref and reads\n")
        val matched = GeneReadsCoupler(
                ref, bam,
                overlap = parm.getOptionValue("min-ref-read", "90.0").toDouble(),
                foldChange = parm.getOptionValue("fold-change", "1.5").toDouble(),
                distanceError = parm.getOptionValue("min-AS-length", "3").toInt()
        )


        println()
        matched.saveTo(File(outDir, "join_stats.tsv").toString())
        matched.saveTemplate(File(outDir, "templates.tsv").toString())

        logger.info("Start to identify splice events\n")
        SJFinder(
                matched,
                distance = parm.getOptionValue("min-AS-length", "3").toInt(),
                overlap = parm.getOptionValue("min-exon-intron", "90.0").toDouble()
        ).saveTo(File(outDir, "final.tsv").toString())


    } catch (err: ParseException) {
        println("Parameter error: $err ")
        help.printHelp("usage message", options)

    } catch (err: IOException) {
        println(err)
    } catch (err: Exception) {
        logger.error(err.message)

        for (i in err.stackTrace) {
            logger.error(i)
        }
    }


}