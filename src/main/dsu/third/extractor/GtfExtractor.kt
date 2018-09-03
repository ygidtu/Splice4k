package dsu.third.extractor

import java.io.File
import java.io.IOException
import java.util.Scanner

import kotlin.system.exitProcess
import org.apache.log4j.Logger

import dsu.carrier.Genes
import dsu.carrier.Exons
import dsu.progressbar.ProgressBar


/**
 * @since 2018.06.14
 * @version 20180903
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

        val msgWeNeed = arrayOf("gene_name", "gene_id", "transcript_name", "transcript_id")

        for (i in 0..(info.size - 1) step 2) {
            if (info[i] !in msgWeNeed) {
                continue
            }
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
        val exons: MutableMap<String, List<Exons>> = mutableMapOf()

        var reader = Scanner(System.`in`)
        val pb = ProgressBar(message = "Reading Gtf")
        try {
            reader = Scanner(File(this.gtf))

            while (reader.hasNext()) {


                pb.step()

                val line = reader.nextLine()

                if (line.startsWith("#")) continue // skip annotations

                val lines = line.split("\\s+".toRegex())
                val tmpGene = Genes(
                        chromosome = lines[0],
                        start = lines[3].toInt(),
                        end = lines[4].toInt(),
                        strand = lines[6].toCharArray()[0],
                        information = this.extractTagInformation(lines.subList(8, lines.size))
                )

                if (lines[2] == "transcript") {
                    transcripts.add(tmpGene)
                } else if (lines[2] == "exon") {
                    val tmp = mutableListOf(Exons(tmpGene.start, tmpGene.end))

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
