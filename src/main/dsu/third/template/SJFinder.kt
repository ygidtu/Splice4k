package dsu.third.template

import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.io.PrintWriter
import java.util.Objects
import dsu.carrier.Exons
import dsu.carrier.SpliceJunction
import dsu.progressbar.ProgressBar



/**
 * @since 2018.06.21
 * @version 20180903
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */


class SJFinder(
        private val pair: GeneReadsCoupler,
        private val distance: Int = 3,
        private val overlap: Double = 90.0
) {
    private val logger = Logger.getLogger(SJFinder::class.java.toString())
    private val results = mutableListOf<SpliceJunction>()

    init {
        this.identifySJ()
    }

    /**
     * 识别各种可变剪接类型
     */
    private fun identifySJ() {


        // 这个gap就是为了控制输出一个合适的进度条的
        val pb = ProgressBar(this.pair.templates.size.toLong(), "Splice events identifying at")
        for (pair in this.pair.templates) {

            pb.step()

            val splice = SpliceJunction(pair.reference)

            this.compareSites(
                    reference = pair.geneExons,
                    reads = pair.template.exons,
                    splice = splice
            )

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
     * @param reference 基因外显子构成的列表
     * @param reads Reads外显子构成的列表
     * @param splice 收集可变剪接类型的类
     * @return SpliceJuntion 某基因上所有可变剪接的类型
     */
    private fun compareSites(
            reference: List<Exons>,
            reads: List<Exons>,
            splice: SpliceJunction
    ) {
        var i = 0; var j = 0
        var skipped = false

        val spliced = mutableSetOf<Int>()

        while (i < reference.size && j < reads.size) {
            val tmpRef = reference[i]
            val tmpRead = reads[j]

            when {
                tmpRef.isUpStream(tmpRead) -> {  // 基因在上游
                    /*
                     如果基因的外显子，并没有与任何read的外显子有重合，
                     那么这个exon基本就是exon_skipping了

                     2018.07.03 修正
                     如果基因已经在上游了，
                     要么他是最开始，
                     要么就已经与合适范围内的所有reads比对过一遍了
                     这种情况下，它还没有任何重合的reads，就必定是skipping了
                      */
                    if (Objects.hash(tmpRef) !in spliced) {
                        splice.addEvent("exon_skipping", tmpRef.start, tmpRef.end)
                        spliced.add(Objects.hash(tmpRef))
                    }


                    /*
                    exon inclusion
                    reads的exon一定在注释的intron范围内
                    因此，
                     */
                    if (
                            i < reference.size - 1 &&
                            tmpRead.start > tmpRef.end &&
                            tmpRead.end < reference[i+1].start
                    ) {
                        splice.addEvent("exon_inclusion", tmpRead.start, tmpRead.end)
                        spliced.add(Objects.hash(tmpRef))
                    }

                    i++

                }

                tmpRef.isDownStream(tmpRead) -> {    // 基因在下游
                    j++
                }

                else -> {   // 有重合

                    spliced.add(Objects.hash(tmpRef))

                    if (  // intron retention   require 90% coverage
                            i < reference.size -1 &&
                            tmpRef.end < reference[i + 1].start &&
                            tmpRead.overlapPercent(
                                    Exons(tmpRef.end, reference[i+1].start)
                            ) > this.overlap
                    ) {
                        splice.addEvent("intron_retention", tmpRef.end, reference[i+1].start)

                        /*
                        由于intron retention
                        reads的该外显子可能横跨多个基因的外显子，
                        因此，通过一个循环，找出其横跨的最后一个基因外显子的index是多少

                        如果是横跨到基因的最后一个外显子了，那么index肯定为最后一个

                        如果不是，那么index需要减一，才是他横跨的最后一个外显子的index

                        再次基础上调整临时基因的范围，去判断是否存在donor/acceptor
                         */
                        var k = i
                        while (k < reference.size && tmpRead.end > reference[k].start) {
                            k++
                        }
                        k--
                        skipped = true
                        i = k + 1
                        tmpRef.end = reference[k].end
                    }



                    if (  // intron in exon
                        j < reads.size - 1 &&
                                tmpRead.end > tmpRef.start &&
                                reads[j + 1].start < tmpRef.end
                    ) {
                        splice.addEvent("intron_in_exon", tmpRead.end, reads[j + 1].start)

                        /*
                        intron in exon就表明，
                        基因的外显子会横跨多个reads的外显子

                        因此，同上，通过循环找出横跨的最后一个reads外显子的index
                         */
                        var k = j
                        while (k < reads.size && tmpRef.end > reads[k].start) {
                            k++
                        }
                        k--
                        skipped = true
                        j = k + 1
                        tmpRead.end = reads[k].end
                    }

                    if (j != 0 && j != reads.size - 1) {
                        // 在不同情形下，确定了不同的基因和reads外显子范围，用来比对donor/acceptor
                        when {
                            !this.isSameRegion(tmpRef.start, tmpRead.start) &&
                                    !this.isSameRegion(tmpRef.end, tmpRead.end) ->
                                splice.addEvent(
                                        "donor/acceptor",
                                        tmpRead.start,
                                        tmpRead.end
                                )

                            !this.isSameRegion(tmpRef.start, tmpRead.start) &&
                                    this.isSameRegion(tmpRef.end, tmpRead.end) ->
                                splice.addEvent(
                                        "donor",
                                        kotlin.math.min(tmpRef.start, tmpRead.start),
                                        kotlin.math.max(tmpRef.start, tmpRead.start)
                                )

                            this.isSameRegion(tmpRef.start, tmpRead.start) &&
                                    !this.isSameRegion(tmpRef.end, tmpRead.end) ->
                                splice.addEvent(
                                        "acceptor",
                                        kotlin.math.min(tmpRef.end, tmpRead.end),
                                        kotlin.math.max(tmpRef.end, tmpRead.end)
                                )
                        }
                    }

                    if (!skipped) {
                        j++
                    } else {
                        skipped = false
                    }
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
