package main.template

import java.io.File
import java.io.PrintWriter
import java.io.IOException
import org.apache.log4j.Logger

import main.carrier.*
import main.extractor.*

/**
 * @author zhangyiming
 * @since 2018.06.20
 * @version 0.1
 * 将基因与reads匹配到一起
 */


class GeneReadsCoupler(private val Gene: Extractor, private val Reads: Extractor, private val silent: Boolean = false) {
    private val logger = Logger.getLogger(GeneReadsCoupler::class.java)

//    val matchedGeneRead = this.matchGeneReads()
    val matchedGeneRead = mutableListOf<GeneRead>()
    val novelReads = mutableListOf<Genes>()

    init {
        this.matchGeneReads()
    }

    /**
     * 将基因与reads匹配到一起
     */
    private fun matchGeneReads() {
        var logged = -1
        val tmpMatched = mutableMapOf<Genes, MutableList<GeneRead>>()

        var firstOverlap = true
        var readIndex = 0
        var tmpGene = this.Gene.next()
        var tmpRead = this.Reads.next()

        // 统计所有的配对信息
        while (tmpGene != null && tmpRead != null) {
            if (this.Gene.index % 10000 == 0) {
                if (logged != this.Gene.index) {
                    logger.info("Gene Reads matching at ${this.Gene.index}/${this.Gene.totalLine}")
                    logged = this.Gene.index
                }

            }

            when {
                // 基因在read上游，下一个基因
                tmpGene.isUpStream(tmpRead) -> {
                    tmpGene = this.Gene.next()

                    // rollback read index
                    this.Reads.index = readIndex
                    tmpRead = this.Reads.get(readIndex)
                    firstOverlap = true
                }
                // 基因在read下游，读一个read
                tmpGene.isDownStream(tmpRead) -> {
                    if (!tmpMatched.containsKey(tmpRead) && tmpRead.exons.size == 1) {
                        this.novelReads.add(tmpRead)
                    }

                    tmpRead = this.Reads.next()
                }

                else -> {
                    // log read index
                    if (firstOverlap) {
                        readIndex = this.Reads.index
                        firstOverlap = false
                    }
                    val tmpGeneRead = GeneRead(tmpGene, tmpRead)

                    // 判断是否临时的匹配中是否含有该条read了
                    val tmpList = mutableListOf(tmpGeneRead)

                    if (tmpMatched.containsKey(tmpRead)) {
                        tmpList.addAll(tmpMatched[tmpRead]!!)
                    }

                    tmpMatched[tmpRead] = tmpList

                    tmpRead = this.Reads.next()

                }
            }
        }


        for (v in tmpMatched.values) {
            when {
                v.size > 1 -> {
//                    results.addAll(v)
                    val tmpV = v.sortedBy { it.overlap.dec() }

                    if (tmpV[0].overlap / tmpV[1].overlap.toDouble() > 1.5) {
                        this.matchedGeneRead.add(tmpV[0])
                    }
                }

                v.size == 1 -> this.matchedGeneRead.add(v[0])
            }
        }

//        return results

    }

    /**
     * 保存基因与Reads匹配的样本和novel的read到文件
     * @param outfile 输出文件路径
     */
    fun saveTo(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for (i in this.matchedGeneRead) {
                writer.println(i)
            }

            for (i in this.novelReads) {
                writer.println("None|$i")
            }

        } catch (err: IOException) {
            logger.error(err.message)
            for (i in err.stackTrace) {
                logger.error(i)
            }
        } finally {
            writer.close()
        }
    }
}

/*
fun main(args: Array<String>) {
    val gene = GffExtractor("/home/zhang/genome/Homo_sapiens.GRCh38.91.gff3")
    val reads = BamExtractor("/home/zhang/splicehunter_test/test.bam", silent = true)

    val test = GeneReadsCoupler(gene, reads)

    test.saveTo("/home/zhang/splicehunter_test/gene_read1.txt")


}
*/