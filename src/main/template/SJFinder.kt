package main.template

import org.apache.log4j.Logger

import java.io.File
import java.io.IOException
import java.io.PrintWriter

import java.util.Objects

import kotlin.math.*

import main.carrier.SpliceJunction
import main.carrier.isDownStream
import main.carrier.isUpStream
import main.carrier.overlapPercent

import main.extractor.*



/**
 * @since 2018.06.21
 * @version 0.1
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */


class SJFinder(
        private val pair: GeneReadsCoupler,
        private val distance: Int = 3,
        private val overlap: Double = 90.0
) {
    private val logger = Logger.getLogger(SJFinder::class.java)
    private val results = mutableListOf<SpliceJunction>()


    init {
        this.identifySJ()
    }

    /**
     * 识别各种可变剪接类型
     */
    private fun identifySJ() {


        // 这个gap就是为了控制输出一个合适的进度条的
        val gap = this.pair.templates.size.toString().length - 2

        for ((index, pair) in this.pair.templates.withIndex()) {
            if (index % 10.0.pow(gap.toDouble()).toInt() == 0) {
                logger.info("Splice events identifying at $index/${this.pair.templates.size}")
            }

            val splice = SpliceJunction(pair.gene)
            this.compareSites(pair.geneExons, pair.template.exons, splice)
            this.results.add(splice)
        }
    }

    /**
     * 检查两个点是否在特定范围内
     *
     */
    private fun isSameRegion(first: Int, second: Int): Boolean {
        return kotlin.math.abs(first - second) <= this.distance
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
        var firstOverlap = true
        var logged = 0

        val spliced = mutableSetOf<Int>()

        while (i < gene.size && j < reads.size) {
            val tmpGene = gene[i]
            val tmpRead = reads[j]

            when {
                isUpStream(tmpGene, tmpRead) -> {  // 基因在上游
                    /*
                     如果基因的外显子，并没有与任何read的外显子有重合，
                     那么这个exon基本就是exon_skipping了

                     2018.07.03 修正
                     如果基因已经在上游了，
                     要么他是最开始，
                     要么就已经与合适范围内的所有reads比对过一遍了
                     这种情况下，它还没有任何重合的reads，就必定是skipping了
                      */
                    if (Objects.hash(tmpGene) !in spliced) {
                        splice.addEvent("exon_skipping", tmpGene[0], tmpGene[1])
                        spliced.add(Objects.hash(tmpGene))
                    }


                    /*
                    exon inclusion
                    reads的exon一定在注释的intron范围内
                    因此，
                     */
                    if (
                            i < gene.size - 1 &&
                            tmpRead[0] > tmpGene[1] &&
                            tmpRead[1] < gene[i+1][0]
                    ) {
                        splice.addEvent("exon_inclusion", tmpRead[0], tmpRead[1])
                        spliced.add(Objects.hash(tmpGene))
                    }

                    i++

                    if (!firstOverlap) {
                        j = logged
                        firstOverlap = true
                    }
                }

                isDownStream(tmpGene, tmpRead) -> {    // 基因在下游
                    j++
                }

                else -> {   // 有重合

                    if (firstOverlap) {
                        logged = j
                        firstOverlap = false
                    }
                    spliced.add(Objects.hash(tmpGene))

                    if (  // intron retention   require 90% coverage
                            i < gene.size -1 &&
                            overlapPercent(
                                tmpRead,
                                arrayOf(tmpGene[1], gene[i+1][0])
                            ) > this.overlap
                    ) {
                        splice.addEvent("intron_retention", tmpRead[1], reads[j+1][0])

                        /*
                        由于intron retention
                        reads的该外显子可能横跨多个基因的外显子，
                        因此，通过一个循环，找出其横跨的最后一个基因外显子的index是多少

                        如果是横跨到基因的最后一个外显子了，那么index肯定为最后一个

                        如果不是，那么index需要减一，才是他横跨的最后一个外显子的index

                        再次基础上调整临时基因的范围，去判断是否存在donor/acceptor
                         */
                        var k = i
                        while (k < gene.size && tmpRead[1] > gene[k][0]) {
                            k++
                        }
                        k--
                        tmpGene[1] = gene[k][1]

                        /*
                        同理，检查这个reads是否横跨前边的基因外显子，
                        然后构成一个完成的基因范围
                         */
                        k = i
                        while ( k >= 0 && tmpGene[1] > tmpRead[0]) {
                            k--
                        }
                        if (k < 0) k++
                        tmpGene[0] = gene[k][0]
                    }

                    if (  // intron in exon
                        j < reads.size - 1 &&
                                tmpRead[1] > tmpGene[0] &&
                                reads[j + 1][0] < tmpGene[1]
                    ) {
                        splice.addEvent("intron_in_exon", tmpRead[1], reads[j + 1][0])

                        /*
                        intron in exon就表明，
                        基因的外显子会横跨多个reads的外显子

                        因此，同上，通过循环找出横跨的最后一个reads外显子的index
                         */
                        var k = j
                        while (k < reads.size && tmpGene[1] > reads[k][0]) {
                            k++
                        }
                        k--
                        tmpRead[1] = reads[k][1]

                        /*
                        同上
                         */
                        k = j
                        while (k >= 0 && tmpGene[0] > reads[k][1]) {
                            k--
                        }
                        if (k < 0) k++
                        tmpRead[0] = reads[k][0]
                    }

                    // 在不同情形下，确定了不同的基因和reads外显子范围，用来比对donor/acceptor
                    when {
                        !this.isSameRegion(tmpGene[0], tmpRead[0]) &&
                                !this.isSameRegion(tmpGene[1], tmpRead[1]) ->
                            splice.addEvent(
                                    "donor/acceptor",
                                    tmpRead[0],
                                    tmpRead[1]
                            )

                        !this.isSameRegion(tmpGene[0], tmpRead[0]) &&
                                this.isSameRegion(tmpGene[1], tmpRead[1]) ->
                            splice.addEvent(
                                    "donor",
                                    kotlin.math.min(tmpGene[0], tmpRead[0]),
                                    kotlin.math.max(tmpGene[0], tmpRead[0])
                            )

                        this.isSameRegion(tmpGene[0], tmpRead[0]) &&
                                !this.isSameRegion(tmpGene[1], tmpRead[1]) ->
                            splice.addEvent(
                                    "acceptor",
                                    kotlin.math.min(tmpGene[1], tmpRead[1]),
                                    kotlin.math.max(tmpGene[1], tmpRead[1])
                            )
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

   test.saveTo("/home/zhang/splicehunter_test/stat1.txt")

   println(Regex("\\.gff3?$").containsMatchIn("/home/zhang/genome/Homo_sapiens.GRCh38.91.gff3"))

}
*/