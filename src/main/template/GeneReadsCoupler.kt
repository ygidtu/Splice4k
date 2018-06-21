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


class GeneReadsCoupler(private val Gene: Extractor, private val Reads: Extractor) {
    private val logger = Logger.getLogger(GeneReadsCoupler::class.java)

    private val matchedGeneRead = this.matchGeneReads()
    /**
     * 将基因与reads匹配到一起
     */
    private fun matchGeneReads(): List<GeneRead> {
        val discardReads = mutableSetOf<Genes>()
        val tmpMatched = mutableMapOf<Genes, GeneRead>()

        var firstOverlap = true
        var readIndex = 0
        var tmpGene = this.Gene.next()
        var tmpRead = this.Reads.next()


        while (tmpGene != null && tmpRead != null) {
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
                    // 保证读取的read不是discarded的
                    while (tmpRead != null && discardReads.contains(tmpRead)) {
                        tmpRead = this.Reads.next()
                    }
                }

                else -> {
                    // log read index
                    if (firstOverlap) {
                        readIndex = this.Reads.index
                        firstOverlap = false
                    }
                    val tmpGeneRead = GeneRead(tmpGene, tmpRead)

                    // 判断是否临时的匹配中是否含有该条read了
                    when(tmpMatched.containsKey(tmpRead)) {
                        true -> {  // 有，判断一下重合程度
                            val overlap = arrayOf(
                                    tmpGeneRead.overlap,
                                    tmpMatched[tmpRead]!!.overlap
                            )

                            overlap.sort()

                            when(overlap[1] / overlap[0].toDouble() < 1.5) {
                                true -> {
                                    tmpMatched.remove(tmpRead)
                                    discardReads.add(tmpRead)
                                }
                                false -> {
                                    when(tmpGeneRead.overlap > tmpMatched[tmpRead]!!.overlap) {
                                        true -> tmpMatched[tmpRead] = tmpGeneRead
                                    }
                                }
                            }
                        }
                        false -> {
                            tmpMatched[tmpRead] = tmpGeneRead
                        }
                    }

                    // 保证读取的read不是discarded的
                    while (tmpRead != null && discardReads.contains(tmpRead)) {
                        tmpRead = this.Reads.next()
                    }
                }
            }
        }

        return tmpMatched.values.toList()

    }

    fun saveTo(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for (i in this.matchedGeneRead) {
                writer.println(i)
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


fun main(args: Array<String>) {
    val gene = GffExtractor("/home/zhang/genome/Homo_sapiens.GRCh38.91.gff3")
    val reads = BamExtractor("/home/zhang/splicehunter_test/test.bam")


}
