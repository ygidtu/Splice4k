package com.splice4k.tools

import com.splice4k.base.Sites
import com.splice4k.index.SJIndex
import java.io.File
import java.io.PrintWriter


/**
 * 用于提取输出每一个样本中，所有junctions的psi table，分成两个表，same start和same end
 * @author Zhang Yiming
 * @since 2018.10.23
 * @version 20181023
 */


class PSITable {
    private var total = 0
    val sameStart = mutableMapOf<String, MutableMap<String, String>>()
    val sameEnd = mutableMapOf<String, MutableMap<String, String>>()
    val samples = mutableListOf<String>()

    /**
     * 计算所有情况下的PSI
     * @param values 某个junctionGraph中记录的所有位点信息
     * @param collection 用于收集所有PSI值的map，即this.sameStart this.sameEnd
     * @param chromosome 染色体
     * @param strand 正负链
     */
    private fun calculateAllPSI(
            values: Iterable<Sites>,
            collection: MutableMap<String, MutableMap<String, String>>,
            chromosome: String,
            strand: Char
    ) {

        for ( site in values ) {
            for (i in site.getSites() ) {
                val key = "$chromosome:${site.node}-${i.site}$strand"

                val listOfPSI = collection[key] ?: mutableMapOf()

                listOfPSI[this.samples.last()] = site.getPsi(i.site).toString()

                collection[key] = listOfPSI
            }
        }

    }


    /**
     * 添加某次计算过程中收集到的所有junctions等信息
     * @param index
     */
    fun addJunctionGraph(index: SJIndex) {
        this.total++

        this.samples.add( index.infile.name )

        for ( graph in index.data.values ) {
            this.calculateAllPSI(
                    values = graph.starts.values,
                    collection = this.sameStart,
                    chromosome = graph.chromosome,
                    strand = graph.strand
            )

            this.calculateAllPSI(
                    values = graph.ends.values,
                    collection = this.sameEnd,
                    chromosome = graph.chromosome,
                    strand = graph.strand
            )
        }
    }


    /**
     * 将所有的psi写出到文件
     * @param prefix 输出文件的前缀名
     */
    fun writeTo(  prefix: File ) {

        fun write(
                output: String,
                collection: MutableMap<String, MutableMap<String, String>>
        ) {
            val writer = PrintWriter(File(output))

            writer.println("#junctions\t${this.samples.joinToString(separator = "\t")}")

            for ( (junction, value) in collection) {
                writer.println("$junction\t${this.samples
                        .asSequence()
                        .map { value[it] ?: "0" }
                        .toList()
                        .joinToString(separator = "\t")}")
            }

            writer.close()
        }

        write( "${prefix}.sameStart.PSI.tab", this.sameStart )

        write( "${prefix}.sameEnd.PSI.tab", this.sameEnd )
    }
}
