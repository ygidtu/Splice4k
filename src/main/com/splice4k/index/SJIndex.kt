package com.splice4k.index

import org.apache.log4j.Logger
import com.splice4k.base.SpliceGraph
import com.splice4k.progressbar.ProgressBar
import java.io.File
import java.util.*
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since ??
 * @version 20180927
 */


/**
 * 读取提取Splice4k得到的SJ文本格式，或者进一步支持STAR的out.tab格式
 * @param infile 输入文件
 * @param fielter 过滤低丰度的SJ
 * @param star 是否为STAR的输出结果
 */
open class SJIndex(
        val infile: String,
        val filter: Int,
        val star: Boolean = false
) {
    open val logger = Logger.getLogger(SJIndex::class.java)

    val data = mutableMapOf<String, SpliceGraph>()

    init {
        this.getAllSJ()
    }

    open fun getAllSJ() {
        logger.info("Reading from ${this.infile}")

        val pattern = when ( star ) {
            false -> "^([\\w\\.]+):(\\d+)-(\\d+)([+-\\.]?)\t(\\d+)$".toRegex()
            else -> "^([\\w\\.]+)\\s(\\d+)\\s(\\d+)\\s([12])\\s\\d+\\s\\d+\\s+(\\d+).*$".toRegex()
        }

        val reader = Scanner(File(this.infile))

        val pb = ProgressBar(message = "Reading SJ")

        while (reader.hasNext()) {
            val line = reader.nextLine()
            pb.step()
            try{
                var (chromosome, tmpStart, tmpEnd, strand, count) = pattern.find(line)!!.destructured

                strand = when( strand  ) {
                    "" -> "."
                    "1" -> "+"
                    "2" -> "-"
                    else -> strand
                }

                if ( count.toInt() < this.filter ) {
                    continue
                }

                val key = "$chromosome$strand"
                val tmpGraph = when ( this.data.containsKey(key) ) {
                    true -> this.data[key]!!
                    else -> SpliceGraph(chromosome = chromosome, strand = strand.toCharArray().first())
                }

                tmpGraph.addEdge(tmpStart.toInt(), tmpEnd.toInt(), count.toInt(), count.toInt())

                this.data[key] = tmpGraph
            } catch ( e: NullPointerException ) {
                continue
            }
        }
        pb.close()

        if ( this.data.isEmpty() ) {
            this.logger.error("${this.infile} format error")
            exitProcess("format error".toInt())
        }
    }

}