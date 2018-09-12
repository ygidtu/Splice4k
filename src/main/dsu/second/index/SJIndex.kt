package dsu.second.index

import java.util.Scanner
import java.io.File

import dsu.carrier.Exons
import org.apache.log4j.Logger

open class SJIndex(val infile: String) {
    val logger = Logger.getLogger(SJIndex::class.java)

    val data = mutableMapOf<Exons, Int>()
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

                val tmpExon = Exons(
                        chromosome = chromosome,
                        start = tmpStart.toInt(),
                        end = tmpEnd.toInt(),
                        strand = strand.toCharArray()[0]
                )

                this.data[tmpExon] = count.toInt()

                val start = tmpExon.getSiteHash(true)
                val end = tmpExon.getSiteHash(false)
                when (this.sameStart.containsKey(start)) {
                    true -> this.sameStart[start]!!.add(tmpExon)
                    false -> this.sameStart[start] = mutableListOf(tmpExon)
                }

                when (this.sameEnd.containsKey(end)) {
                    true -> this.sameEnd[end]!!.add(tmpExon)
                    false -> this.sameEnd[end] = mutableListOf(tmpExon)
                }
                
            } catch ( e: NullPointerException ) {
                continue
            }
        }
    }
}