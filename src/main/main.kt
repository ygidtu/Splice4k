package main

import java.io.File
import java.io.IOException

import kotlin.system.exitProcess

import org.apache.commons.cli.Options
import org.apache.commons.cli.DefaultParser
import org.apache.commons.cli.HelpFormatter
import org.apache.commons.cli.ParseException

import org.apache.log4j.Logger

import main.extractor.*
import main.template.*
import kotlin.math.log


/**
 * @author zhangyiming
 * @since 2018.06.23
 * @version 0.1
 *
 * command line parameters
 */

fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")

    val options = Options()

    options.addOption("b", "bam",true, "input bam file")
    options.addOption("g",true, "[gff3|gtf]")
    options.addOption("o", "output", true, "output directory")
    options.addOption("h", "help", false, "print help message")
    options.addOption("s", "silent", false, "suppress output message")

    val parser = DefaultParser()
    val help = HelpFormatter()
    try{
        val parameters = parser.parse(options, args)


        var silent = false
        when {
            parameters.hasOption("help") || args.isEmpty() -> {
                help.printHelp("usage message", options)
                exitProcess(0)
            }

            parameters.hasOption("silent") -> silent = true
        }


        val outDir = File(parameters.getOptionValue("output"))

        if (!outDir.exists()) outDir.parentFile.mkdirs()

        val bamFile = File(parameters.getOptionValue("bam")).absoluteFile
        val bamTsv = File(outDir, bamFile.name.split(".")[0] + ".tsv")

        val bam = BamExtractor(bamFile.toString(), silent = silent)
        bam.saveTo(bamTsv.toString())


        val gff : Extractor

        if (parameters.getOptionValue("g").endsWith("gff3")) {
            val gffFile = File(parameters.getOptionValue("g")).absoluteFile
            val gffTsv = File(outDir, gffFile.name.split(".")[0] + ".tsv")

            gff = GffExtractor(gffFile.toString(), silent)
            gff.saveTo(gffTsv.toString())
        } else if (parameters.getOptionValue("g").endsWith("gtf")) {
            val gffFile = File(parameters.getOptionValue("g")).absoluteFile
            val gffTsv = File(outDir, gffFile.name.split(".")[0] + "tsv")

            gff = GffExtractor(gffFile.toString(), silent)
            gff.saveTo(gffTsv.toString())
        } else {
            println("Input should be gtf or gff3 format")
            exitProcess(2)
        }


        val matched = GeneReadsCoupler(gff, bam)

        matched.saveTo(File(outDir, "join_stats.tsv").toString())


        SJFinder(matched).saveTo(File(outDir, "final.tsv").toString())


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