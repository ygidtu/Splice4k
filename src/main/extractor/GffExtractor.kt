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
 * 从gff3格式中提取所需信息
 */


/**
 * @param gff 输入的gff文件的路径
 * @param silent Boolean值，减少信息输出，默认为false
 */
class GffExtractor(
        private val gff: String,
        private val silent: Boolean = false
) : Extractor(silent) {

//    override val logger = Logger.getLogger(GffExtractor::class.java)

    init {
        this.data = gffReader()
        this.totalLine = this.data.size
    }

    /**
     * 提取每一行gtf中gene_id等标签信息
     * @info 含有标签信息的列表
     * @return 字典<信息类别，信息本身>
     */
    private fun extractTagInformation(info: String): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (inf in info.split(";")) {
            val tmp = inf.split("=")

            if (tmp[0] == "Parent") {
                results["ParentType"] = tmp[1].split(":")[0]
                results["Parent"] = tmp[1].split(":")[1]
            } else {
                results[tmp[0]] = tmp[1]
            }
        }

        return results
    }


    /**
     * 读取整个gtf文件，并提取封装其中信息
     * @return 一列排序过后的Genes
     */
    private fun gffReader(): List<Genes> {
        val transcripts: MutableList<Genes> = mutableListOf()
        val exons: MutableMap<String, List<Array<Int>>> = mutableMapOf()

        var reader = Scanner(System.`in`)
        try{
            reader = Scanner(File(this.gff))
            val pattern = Regex(".*transcript:.*")

            var readed = 0
            while (reader.hasNext()) {
                if (readed % 100000 == 0) {
                    this.logger.info("Read $readed lines")
                }
                readed ++

                val line = reader.nextLine()

                if (!pattern.matches(line)) continue   // skip annotations

                val lines = line.split("\t")

                val tmpGene = Genes(
                        chrom = lines[0],
                        start = lines[3].toInt(),
                        end = lines[4].toInt(),
                        strand = lines[6].toCharArray()[0],
                        information = this.extractTagInformation(lines[8])
                )

                // 分别判断是否为转录本和exon
                if (".*ID=transcript:.*".toRegex().matches(line)) {
                    transcripts.add(tmpGene)
                } else if (".*Parent=transcript:.*".toRegex().matches(line) && lines[2] == "exon") {

                    // 外显子收集的这个写法比Python复杂些，但是功能是一样的
                    val tmp = mutableListOf(arrayOf(tmpGene.start, tmpGene.end))

                    if (exons.containsKey(tmpGene.parent)) {
                        tmp.addAll(exons[tmpGene.parent]!!)
                    }
                    exons[tmpGene.parent] = tmp
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

        // 为转录本指定外显子
        for (i in transcripts) {
            if (exons.containsKey(i.transcriptId)) {
                i.exons = exons[i.transcriptId]!!.toMutableList()
            }
        }

        return transcripts
    }
}

/*
fun main(args: Array<String>) {
    val test = GffExtractor("/home/zhang/genome/Homo_sapiens.GRCh38.91.gff3")

    test.saveTo("/home/zhang/genome/test/hsa_gff.txt")
}
*/
