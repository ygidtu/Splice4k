package com.splice4k

import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.core.subcommands
import com.splice4k.cli.*
import org.apache.log4j.Logger
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since ???
 * @version 20181006
 */

const val VERSION = "Splice4k version: 1.2.4 -> 2018.11.08"


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters()
            .subcommands(Extract())
            .subcommands(SGS())
            .subcommands(SMS())
            .subcommands(Iso())

    // help message
    if (args.size <= 1 || "-h" in args || "--help" in args ) {
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

    }
}

