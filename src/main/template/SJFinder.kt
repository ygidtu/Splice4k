package main.template

import org.apache.log4j.Logger

import java.io.File
import java.io.IOException
import java.io.PrintWriter

import kotlin.math.*

import main.carrier.SpliceJunction
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
     * 检查第一个外显子是否在第二个的上游
     * @param first 外显子位点 Array(start, end)
     * @param second 外显子位点 Array(start, end)
     * @return true 有。false 没有
     */
    private fun isUpStream(first: Array<Int>, second: Array<Int>): Boolean {
        return first[1] < second[0]
    }


    /**
     * 检查第一个外显子是否在第二个的下游
     * @param first 外显子位点 Array(start, end)
     * @param second 外显子位点 Array(start, end)
     * @return true 有。false 没有
     */
    private fun isDownStream(first: Array<Int>, second: Array<Int>): Boolean {
        return first[0] > second[1]
    }


    /**
     * 检查两个点是否在特定范围内
     *
     */
    private fun isSameRegion(first: Int, second: Int): Boolean {
        return kotlin.math.abs(first - second) <= this.distance
    }

    /**
     * 计算两个文件之间的重合程度
     * @param first 第一个位点的坐标，主要是reads上的exon
     * @param second 第二个位点的坐标，主要是基因的intron
     * @return double, 两个位点重合的比例，如果<=0，则没有重合
     */
    private fun overlapPercent(first: Array<Int>, second: Array<Int>): Double {
        return (min(first[1], second[1]) - max(first[0], second[0])) /
                (second[1] - second[0]).toDouble()
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
                this.isUpStream(tmpGene, tmpRead) -> {  // 基因在上游
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
                        spliced.add(tmpGene)
                    }

                    i++
                }

                this.isDownStream(tmpGene, tmpRead) -> {    // 基因在下游
                    /*
                     如果基因的外显子，并没有与任何read的外显子有重合，
                     那么这个exon基本就是exon_skipping了
                      */
                    if (!spliced.contains(tmpGene)) {
                        splice.addEvent("exon_skipping", tmpGene[0], tmpGene[1])
                        spliced.add(tmpGene)
                    }

                    j++
                }

                else -> {   // 有重合

                    when {
                        /*
                        intron_retention
                        应当是reads的exon与intron区域有重合
                        重合比例应当达到90%以上
                         */
                        i < gene.size - 1 &&
                            this.overlapPercent(
                                    tmpRead,
                                    arrayOf(tmpGene[1], gene[i+1][0])
                            ) > this.overlap -> {
                            splice.addEvent("intron_retention", tmpRead[1], reads[j+1][0])
                        }

                        /*
                        如果不属于intron_retention
                        那么就需要判断是否属于donor或者acceptor
                        位点之间距离误差在3以内的就认为是同一个位点
                         */
                        !this.isSameRegion(tmpGene[0], tmpRead[0]) && !this.isSameRegion(tmpGene[1], tmpRead[1]) -> {
                            splice.addEvent("donor/acceptor", tmpGene[0], tmpGene[1])
                        }

                        this.isSameRegion(tmpGene[0], tmpRead[0]) && !this.isSameRegion(tmpGene[1], tmpRead[1]) -> {
                            splice.addEvent("doner", tmpRead[0], tmpRead[1])
                        }

                        !this.isSameRegion(tmpGene[0], tmpRead[0]) && this.isSameRegion(tmpRead[1], tmpRead[1]) -> {
                            splice.addEvent("acceptor", tmpRead[0], tmpRead[1])
                        }

                    }
                    spliced.add(tmpGene)
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