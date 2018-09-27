package splice4k.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import splice4k.index.BamIndex

/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20180927
 */


class Extract: CliktCommand(help = "Extract junctions from Bam/Sam files") {
    private val input by argument(help = "Path to input Bam/Sam file").file(exists = true)
    private val output by argument(help = "Path to output file").file()
    private val junctionsFilter by option(
            "-c",
            "--count",
            help = "Filter low abundance junctions [default: 10]"
    ).int().default(10).validate {
        it >= 0
    }


    override fun run() {
        BamIndex(
                infile = this.input.toString(),
                filter = this.junctionsFilter
        ).writeTo(output = this.output)
    }
}