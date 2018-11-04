package com.splice4k.tools

import com.splice4k.index.SJIndex
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter


/**
 * 用于提取输出每一个样本中，所有junctions的psi table，分成两个表，same start和same end
 * @author Zhang Yiming
 * @since 2018.10.23
 * @version 20181023
 */


class CountTable {
    private val logger = Logger.getLogger(CountTable::class.java)

    /**
     * as function name says
     * format all the junctions counts across all samples into matrix
     * @return Map<Sample, Map<Junction, Count>
     */
    private fun formatJunctionGraphToMap( data: List<SJIndex> ): Map<String, Map<String, Int>> {
        val res = mutableMapOf<String, MutableMap<String, Int>>()
        data.forEach {
            for ( (key, count) in it.data ) {

                val tmpMap = res[it.infile.name] ?: mutableMapOf()

                tmpMap[key] = count
                res[it.infile.name] = tmpMap
            }
        }

        return res
    }

    /**
     * save count table to file
     * @param prefix prefix of output file
     */
    fun writeTo( prefix: File, data: List<SJIndex>) {
        this.logger.info("Extract Junction's Counts")

        fun write(
                data: Map<String, Map<String, Int>>,
                output: String
        ) {
            val writer = PrintWriter(File(output))
            val samples = data.keys.toList()
            val junctions = mutableSetOf<String>()

            data.values.forEach {
                junctions.addAll(it.keys)
            }

            writer.println("#junctions\t${samples.joinToString(separator = "\t")}")


            for ( row in junctions ) {
                var count = ""
                for ( col in samples ) {
                    if ( count != "" ) {
                        count += "\t"
                    }
                    count += data[col]!![row]?.toString() ?: "0"
                }
                val genomicLoci = row.split("\t")
                val genomicLociFormatted = "${genomicLoci[0]}:${genomicLoci[1]}-${genomicLoci[2]}${genomicLoci[3]}"
                writer.println("$genomicLociFormatted\t$count")
            }

            writer.close()
        }

        write(
                data = this.formatJunctionGraphToMap(data),
                output = "$prefix.junction_counts.tab"
        )
    }

    /**
     * collect the junctions that not overall counts lower than threshold
     */
    fun filter(data: List<SJIndex>, threshold: Int): HashSet<String> {
        this.logger.info("Filtering by junctions across samples")

        val res = hashSetOf<String>()

        val junctions = mutableSetOf<String>()

        data.forEach {
            junctions.addAll(it.data.keys)
        }

        this.logger.info("Before: Total junctions: ${junctions.size}")

        for ( junction in junctions ) {
            var totalCount = 0
            data.forEach {
                totalCount += (it.data[junction] ?: 0)
            }

            if ( totalCount < threshold ) {
                res.add(junction)
            }
        }

        this.logger.info("After: Total junctions: ${junctions.size - res.size}")
        return res
    }
}