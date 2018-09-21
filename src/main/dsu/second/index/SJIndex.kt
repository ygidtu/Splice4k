package dsu.second.index

import java.util.Scanner
import java.io.File
import dsu.carrier.SpliceGraph
import org.apache.log4j.Logger
import dsu.progressbar.ProgressBar

open class SJIndex(val infile: String) {
    val logger = Logger.getLogger(SJIndex::class.java)

    val data = mutableMapOf<String, SpliceGraph>()

    init {
        this.getAllSJ()
    }

    open fun getAllSJ() {
        logger.info("Reading Splice Junctions")
        val pattern = "^([\\w\\.]+):(\\d+)-(\\d+)([+-\\.]?)\t(\\d+)$".toRegex()

        val reader = Scanner(File(this.infile))

        val pb = ProgressBar(message = "Reading SJ")

        while (reader.hasNext()) {
            val line = reader.nextLine()
            pb.step()
            try{
                var (chromosome, tmpStart, tmpEnd, strand, count) = pattern.find(line)!!.destructured

                if ( strand == "" ) {
                    strand = "."
                }

                val key = "$chromosome$strand"

                val tmpGraph = when ( this.data.containsKey(key) ) {
                    true -> this.data[key]!!
                    else -> SpliceGraph(chromosome = chromosome, strand = strand.toCharArray().first())
                }

                tmpGraph.addEdge(tmpStart.toInt(), tmpEnd.toInt(), count.toInt(), count.toInt())

                this.data["$chromosome$strand"] = tmpGraph
            } catch ( e: NullPointerException ) {
                continue
            }
        }
        pb.close()
    }
}