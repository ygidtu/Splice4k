package main.template

import main.carrier.SpliceJunction
import main.extractor.BamExtractor
import main.extractor.GffExtractor
import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.io.PrintWriter

/**
 * @since 2018.06.21
 * @version 0.1
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */


class SJFinder(private val pair: GeneReadsCoupler) {
    private val logger = Logger.getLogger(GeneReadsCoupler::class.java)
    val results = mutableListOf<SpliceJunction>()


    init {
        this.identifySJ()
    }

    /**
     * 识别各种可变剪接类型
     */
    fun identifySJ() {
        for (pair in this.pair.matchedGeneRead) {
            val splice = SpliceJunction(pair.gene)
            this.compareSites(pair.gene.exons, pair.reads.exons, splice)
            this.results.add(splice)
        }
    }


    /**
     * 比较位点
     * @param gene 基因外显子构成的列表
     * @param reads Reads外显子构成的列表
     * @param splice 收集可变剪接类型的类
     * @return SpliceJuntion 某基因上所有可变剪接的类型
     */
    private fun compareSites(
            gene: List<Array<Int>>,
            reads: List<Array<Int>>,
            splice: SpliceJunction
    ) {
        var i = 0; var j = 0

        val spliced = mutableSetOf<Array<Int>>()

        while (i < gene.size && j < reads.size) {
            val tmpGene = gene[i]
            val tmpRead = reads[j]

            when {
                tmpGene[1] < tmpRead[0] -> {  // 基因在上游
                    if (
                            i < gene.size - 1 &&
                            tmpRead[0] > tmpGene[1] &&
                            tmpRead[1] < gene[i+1][0]
                    ) {
                        splice.addEvent("exon_inclusion", tmpRead[0], tmpRead[1])
                        spliced.add(tmpGene)
                    }

                    i++
                }
                tmpGene[0] > tmpRead[1] -> {    // 基因在下游
                    if (!spliced.contains(tmpGene)) {
                        splice.addEvent("exon_skipping", tmpGene[0], tmpGene[1])
                        spliced.add(tmpGene)
                    }

                    j++
                }

                else -> {   // 有重合

                    when {
                        tmpGene[0] == tmpRead[0] && tmpGene[1] != tmpRead[1] -> {
                            splice.addEvent("doner", tmpRead[0], tmpRead[1])
                            spliced.add(tmpGene)
                        }

                        tmpGene[0] != tmpRead[0] && tmpRead[1] == tmpRead[1] -> {
                            splice.addEvent("acceptor", tmpRead[0], tmpRead[1])
                            spliced.add(tmpGene)
                        }

                        tmpGene[0] != tmpRead[0] && tmpGene[1] != tmpRead[1] -> {
                            splice.addEvent("donor/acceptor", tmpGene[0], tmpGene[1])
                            spliced.add(tmpGene)
                        }

                        j < reads.size - 1 &&
                                tmpGene[0] < tmpRead[1] &&
                                tmpGene[1] > reads[j+1][0] -> {
                            splice.addEvent("intron_retention", tmpRead[1], reads[j+1][0])
                            spliced.add(tmpGene)
                        }
                    }

                    j++

                }

            }

        }


    }


    /**
     * 保存至文件
     * @param outfile 输出文件的路径
     * @return
     */
    fun saveTo(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for (i in this.results) {
                var line = i.next()

                while (line != null) {
                    writer.println(line)
                    line = i.next()
                }
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

    val test = SJFinder(GeneReadsCoupler(gene, reads))

    test.saveTo("/home/zhang/splicehunter_test/stat.txt")
}
*/