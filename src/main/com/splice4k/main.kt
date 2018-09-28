package com.splice4k

import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.core.subcommands
import org.apache.log4j.Logger
import com.splice4k.cli.Extract
import com.splice4k.cli.Parameters
import com.splice4k.cli.SGS
import com.splice4k.cli.SMRT
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since ???
 * @version 20180927
 */


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters().subcommands(Extract()).subcommands(SGS()).subcommands(SMRT())

    // help message
    if (args.size <= 1) {
        val help = when {
            args.isEmpty() -> cmd.getFormattedHelp()
            args[0].toLowerCase() == "smrt" -> SMRT().getFormattedHelp()
            args[0].toLowerCase() == "extract" -> Extract().getFormattedHelp()
            args[0].toLowerCase() == "sgs" -> SGS().getFormattedHelp()
            args[0].toLowerCase() in arrayOf("-v", "--version") -> {
                println("Splice4k version: 20180927")
                exitProcess(0)
            }
            else -> cmd.getFormattedHelp()
        }
        println(help)
        exitProcess(0)
    }

    // help message
    if ( "-h" in args || "--help" in args ) {
        val help = when (args[0].toLowerCase()) {
            "smrt" -> SMRT().getFormattedHelp()
            "extract" -> Extract().getFormattedHelp()
            "sgs" -> SGS().getFormattedHelp()
            else -> cmd.getFormattedHelp()
        }

        println(help)
        exitProcess(0)
    }


    try{
        cmd.parse(args)


        logger.info("Done")

    } catch (e: UsageError) {

        println(e.localizedMessage)

        println("=".repeat(40))
        when ( args[0].toLowerCase() ) {
            "smrt" -> println(SMRT().getFormattedUsage())
            "sgs" -> println(SGS().getFormattedUsage())
            "extract" -> println(Extract().getFormattedUsage())
            else -> println(cmd.getFormattedUsage())
        }

    } catch (e: Exception) {

        logger.error(e.localizedMessage)

        e.stackTrace.forEach { logger.error(it) }

    }
}
