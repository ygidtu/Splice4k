package dsu.second.index

import java.util.Scanner
import java.io.File

import dsu.carrier.Exons
import dsu.carrier.SpliceGraph
import org.apache.log4j.Logger

open class SJIndex(val infile: String) {
    val logger = Logger.getLogger(SJIndex::class.java)

    val data = mutableMapOf<String, SpliceGraph>()
    val sameStart = mutableMapOf<Int, MutableList<Exons>>()
    val sameEnd = mutableMapOf<Int, MutableList<Exons>>()

    init {
        this.getAllSJ()
    }

    open fun getAllSJ() {
        logger.info("Reading Splice Junctions")
        val pattern = "^([\\w\\.]+):(\\d+)-(\\d+)([+-\\.])\t(\\d+)$".toRegex()

        val reader = Scanner(File(this.infile))

        while (reader.hasNext()) {
            val line = reader.nextLine()

            try{
                val (chromosome, tmpStart, tmpEnd, strand, count) = pattern.find(line)!!.destructured

                val key = "$chromosome$strand"

                val tmpGraph = when ( this.data.containsKey(key) ) {
                    true -> this.data[key]!!
                    else -> SpliceGraph(chromosome=chromosome, strand = strand.toCharArray().first())
                }

                tmpGraph.addEdge(tmpStart.toInt(), tmpEnd.toInt(), count.toInt(), count.toInt())

            } catch ( e: NullPointerException ) {
                continue
            }
        }
    }
}