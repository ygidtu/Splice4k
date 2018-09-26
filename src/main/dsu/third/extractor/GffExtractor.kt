package dsu.third.extractor

import dsu.carrier.Exons
import dsu.carrier.Genes
import dsu.progressbar.ProgressBar
import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.util.*
import kotlin.system.exitProcess


/**
 * @since 2018.06.14
 * @version 20180926
 * @author Zhang yiming
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

    private val logger = Logger.getLogger(GffExtractor::class.java)

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

            if ( ":" in tmp[1] ) {
                results[tmp[0]] = tmp[1].split(":")[1]
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
        val exons: MutableMap<String, MutableList<Int>> = mutableMapOf()
        val geneTranscript = mutableMapOf<String, String>()

        val pb = ProgressBar(message = "Reading Gff")
        var reader = Scanner(System.`in`)
        try{
            reader = Scanner(File(this.gff))
            while (reader.hasNext()) {
                pb.step()
                val line = reader.nextLine()

                val lines = line.split("\\s+".toRegex())

                if ( line.startsWith("#") ) {
                    continue
                }

                val info = this.extractTagInformation(lines[8])
                val tmpGene = Genes(
                        chromosome = lines[0],
                        start = lines[3].toInt(),
                        end = lines[4].toInt(),
                        strand = lines[6].toCharArray()[0],
                        information = info
                )

                // 分别判断是否为转录本和exon
                if ( "transcript_id" in info.keys ) {
                    transcripts.add(tmpGene)
                    geneTranscript[info["ID"]!!] = info["Parent"]!!

                } else if ( lines[2] == "exon" ) {
                    // 外显子收集的这个写法比Python复杂些，但是功能是一样的
                    val tmp = mutableListOf(tmpGene.start, tmpGene.end)

                    if (exons.containsKey(tmpGene.parent)) {
                        tmp.addAll(exons[tmpGene.parent]!!)
                    }
                    exons[tmpGene.parent] = tmp


                    val tmpExon = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            exonId = info["exon_id"]!!
                    )

                    tmpExon.source["transcript"] = info["Parent"]!!
                    tmpExon.source["gene"] = geneTranscript[info["Parent"]]!!

                    val tmpIndex = mutableListOf(tmpExon)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.index.containsKey(key) ) {
                        this.index[key]!!.addAll(tmpIndex)
                    } else {
                        this.index[key] = tmpIndex
                    }
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

        // 为转录本指定外显子
        for (i in transcripts) {
            if (exons.containsKey(i.transcriptId)) {
                i.exons.addAll(exons[i.transcriptId]!!)
            }
        }

        return transcripts
    }
}
