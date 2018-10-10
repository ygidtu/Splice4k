package com.splice4k

import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.core.subcommands
import org.apache.log4j.Logger
import com.splice4k.cli.Extract
import com.splice4k.cli.Parameters
import com.splice4k.cli.SGS
import com.splice4k.cli.SMS
import kotlin.system.exitProcess
import com.splice4k.cli.Iso


/**
 * @author Zhang Yiming
 * @since ???
 * @version 20181006
 */

const val VERSION = "Splice4k version: 1.1.0 -> 2018.10.10"


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters()
            .subcommands(Extract())
            .subcommands(SGS())
            .subcommands(SMS())
            .subcommands(Iso())

    // help message
    if (args.size <= 1) {
        val help = when {
            args.isEmpty() -> cmd.getFormattedHelp()
            args[0].toLowerCase() == "sms" -> SMS().getFormattedHelp()
            args[0].toLowerCase() == "extract" -> Extract().getFormattedHelp()
            args[0].toLowerCase() == "sgs" -> SGS().getFormattedHelp()
            args[0].toLowerCase() == "iso" -> Iso().getFormattedHelp()
            args[0].toLowerCase() in arrayOf("-v", "--version") -> {
                println(VERSION)
                exitProcess(0)
            }
            else -> cmd.getFormattedHelp()
        }
        println(help)
        exitProcess(0)
    }

    // help message
    if ( "-h" in args || "--help" in args ) {
        println(VERSION)
        val help = when (args[0].toLowerCase()) {
            "sms" -> SMS().getFormattedHelp()
            "extract" -> Extract().getFormattedHelp()
            "sgs" -> SGS().getFormattedHelp()
            "iso" -> Iso().getFormattedHelp()
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
        println(VERSION)
        when ( args[0].toLowerCase() ) {
            "sms" -> println(SMS().getFormattedHelp())
            "sgs" -> println(SGS().getFormattedHelp())
            "extract" -> println(Extract().getFormattedHelp())
            "iso" -> println(Iso().getFormattedHelp())
            else -> println(cmd.getFormattedHelp())
        }

    } catch (e: Exception) {

        logger.error(e.localizedMessage)

        e.stackTrace.forEach { logger.error(it) }

        exitProcess(e.hashCode())
    }
}

