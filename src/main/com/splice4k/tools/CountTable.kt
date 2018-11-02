package com.splice4k.tools

import com.splice4k.base.JunctionsGraph
import com.splice4k.index.SJIndex
import java.io.File
import java.io.PrintWriter


/**
 * 用于提取输出每一个样本中，所有junctions的psi table，分成两个表，same start和same end
 * @author Zhang Yiming
 * @since 2018.10.23
 * @version 20181023
 */


class CountTable {
    private val junctionCounts = mutableMapOf<String, Map<String, Int>>()
    private val samples = mutableListOf<String>()

    /**
     * 添加某次计算过程中收集到的所有junctions等信息
     * @param index
     */
    fun addJunctionGraph(index: SJIndex) {
        this.samples.add( index.infile.name )

        this.junctionCounts[this.samples.last()] = this.formatJunctionGraphToMap(index.data)
    }

    /**
     * as function name says
     *
     */
    private fun formatJunctionGraphToMap( data: Map<String, JunctionsGraph> ): Map<String, Int> {
        val res = mutableMapOf<String, Int>()
        for ( (key, graph) in data ) {

            val strand = key.toCharArray().last()

            val chromosome = key.toCharArray().toList()

            for ( i in graph ) {
                for ( j in i ) {
                    val site = "${chromosome.subList(0, chromosome.size - 1).joinToString(separator = "")}:" +
                            j.first +
                            "$strand"

                    res[site] = j.second
                }
            }
        }
        return res
    }

    /**
     * save count table to file
     * @param prefix prefix of output file
     */
    fun writeTo( prefix: File) {

        fun write(
                output: String,
                collection: MutableMap<String, Map<String, Int>>
        ) {
            val writer = PrintWriter(File(output))

            writer.println("#junctions\t${this.samples.joinToString(separator = "\t")}")

            for ( (junction, value) in collection) {
                writer.println("$junction\t${this.samples
                        .asSequence()
                        .map { value[it] ?: 0 }
                        .toList()
                        .joinToString(separator = "\t")}")
            }

            writer.close()
        }

        write( "$prefix.junction_counts.tab", this.junctionCounts )
    }
}