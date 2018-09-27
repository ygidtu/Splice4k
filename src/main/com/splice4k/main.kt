package splice4k

import com.github.ajalt.clikt.core.CliktError
import com.github.ajalt.clikt.core.subcommands
import org.apache.log4j.Logger
import splice4k.cli.Extract
import splice4k.cli.Parameters
import splice4k.cli.SGS
import splice4k.cli.SMRT
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since ???
 * @version 20180927
 */


fun main(args: Array<String>) {
    val logger = Logger.getLogger("main")
    val cmd = Parameters().subcommands(Extract()).subcommands(SGS()).subcommands(SMRT())
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

    try{
        cmd.parse(args)
        logger.info("Done")
    } catch (e: CliktError) {
        logger.error(e.localizedMessage)
        println()
        when ( args[0].toLowerCase() ) {
            "smrt" -> println(SMRT().getFormattedHelp())
            "sgs" -> println(SGS().getFormattedHelp())
            "extract" -> println(Extract().getFormattedHelp())
            else -> println(cmd.getFormattedHelp())
        }
    }
}
