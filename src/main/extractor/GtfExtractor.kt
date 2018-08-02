package main.extractor

import java.io.File
import java.io.IOException
import java.util.Scanner

import kotlin.system.exitProcess
import org.apache.log4j.Logger

import main.carrier.Genes


/**
 * @since 2018.06.14
 * @version 0.1
 * @author zhangyiming
 * 从gtf格式中提取所需信息
 */


/**
 * @param gtf 输入的gtf文件的路径
 * @param silent Boolean值，减少信息输出，默认为false
 */
class GtfExtractor(
        private val gtf: String,
        private val silent: Boolean = false
): Extractor(silent) {

    private val logger = Logger.getLogger(GtfExtractor::class.java)

    init {
        this.data = gtfReader()
        this.totalLine = this.data.size
    }

    /**
     * 提取每一行gtf中gene_id等标签信息
     * @info 含有标签信息的列表
     * @return 字典<信息类别，信息本身>
     */
    private fun extractTagInformation(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (i in 0..(info.size - 1) step 2) {
            results[info[i]] = info[i+1]     // get rid of useless characters
                    .replace("\"", "")
                    .replace(";", "")

        }
        return results
    }


    /**
     * 读取整个gtf文件，并提取封装其中信息
     * @return 字典<gene id, Genes>
     */
    private fun gtfReader(): List<Genes> {
        val transcripts: MutableList<Genes> = mutableListOf()
        val exons: MutableMap<String, List<Array<Int>>> = mutableMapOf()

        var reader = Scanner(System.`in`)

        try {
            reader = Scanner(File(this.gtf))

            var readIn = 0
            var gap = 10
            while (reader.hasNext()) {
                if (readIn % gap == 0) {
                    this.logger.info("Read $readIn lines")
                    if (gap < 10001) gap *= 10
                }
                readIn ++

                val line = reader.nextLine()

                if (line.startsWith("#")) continue // skip annotations

                val lines = line.split("\\s+".toRegex())
                val tmpGene = Genes(
                        chrom=lines[0],
                        start=lines[3].toInt(),
                        end=lines[4].toInt(),
                        strand = lines[6].toCharArray()[0],
                        information = this.extractTagInformation(lines.subList(8, lines.size))
                )

                if (lines[2] == "transcript") {
                    transcripts.add(tmpGene)
                } else if (lines[2] == "exon") {
                    val tmp = mutableListOf(arrayOf(tmpGene.start, tmpGene.end))

                    if (exons.containsKey(tmpGene.geneId)) {
                        tmp.addAll(exons[tmpGene.geneId]!!)
                    }
                    exons[tmpGene.geneId] = tmp
                }
            }
        } catch (err: IOException) {
            this.logger.error(err.message)

            for (i in err.stackTrace) {
                this.logger.error(i)
            }

            exitProcess(1)
        } finally {
            reader.close()
        }

        transcripts.sort()

        for (i in transcripts) {
            if (exons.containsKey(i.geneId)) {
                i.exons = exons[i.geneId]!!.toMutableList()
            }
        }

        return transcripts
    }

}


/*
fun main(args: Array<String>) {
    val test = GtfExtractor("/home/zhang/genome/gencode.v19.annotation.gtf")
    test.saveTo("/home/zhang/genome/test/hsa_gtf.txt")


}
*/
