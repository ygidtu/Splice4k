package main.template

import java.util.Scanner
import java.io.File
import java.io.PrintWriter
import java.io.IOException

import kotlin.math.*

import org.apache.log4j.Logger

import main.carrier.*
import main.extractor.*


/**
 * @author zhangyiming
 * @since 2018.06.20
 * @version 0.1
 * 将基因与reads匹配到一起
 */


class GeneReadsCoupler(
        private val Gene: Extractor,
        private val Reads: Extractor,
        private val overlap: Double = 90.0,
        private val foldChange: Double = 1.5,
        private val distanceError: Int = 3
        ) {
    private val logger = Logger.getLogger(GeneReadsCoupler::class.java)

    val matchedGeneRead = mutableListOf<GeneRead>()
    val novelReads = mutableListOf<GeneRead>()

    lateinit var templates : MutableList<Template>

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
                /*
                 基因在read上游，下一个基因
                 添加了3bp的误差空间，如果距离在3bp内的都算是同一个点了
                 */
                tmpGene.isUpStream(tmpRead, this.distanceError) -> {
                    tmpGene = this.Gene.next()

                    // rollback read index
                    this.Reads.index = readIndex
                    tmpRead = this.Reads.get(readIndex)
                    firstOverlap = true
                }
                // 基因在read下游，读一个read
                tmpGene.isDownStream(tmpRead, this.distanceError) -> {
                    if (!tmpMatched.containsKey(tmpRead) && tmpRead.exons.size == 1) {
                        this.novelReads.add(GeneRead(Genes(), tmpRead))
                    }

                    tmpRead = this.Reads.next()
                }

                else -> {
                    if (tmpGene.strand != tmpRead.strand) {
                        this.novelReads.add(GeneRead(tmpGene, tmpRead))
                    } else {
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
                    }

                    tmpRead = this.Reads.next()
                }
            }
        }

        // 添加重合程度判断，相较于两者中短的区域，覆盖程度要大于90%
        for (v in tmpMatched.values) {
            when {
                v.size > 1 -> {
                    val tmpV = v.sortedBy { it.overlap.dec() }

                    if (
                            tmpV[0].overlap / tmpV[1].overlap.toDouble() > this.foldChange &&
                            tmpV[0].overlapPercent > this.overlap
                    ) {
                        this.matchedGeneRead.add(tmpV[0])
                    }
                }

                v.size == 1 -> {
                    if (v[0].overlapPercent > this.overlap) {
                        this.matchedGeneRead.add(v[0])
                    }
                }
            }
        }

        this.matchedGeneRead.sort()
        this.novelReads.sort()

        this.buildTemplates()
    }


    /**
     * 组装template
     * 有几点要求，
     * 1. 属于同一个基因的reads统统取来做组装
     * 2. 先将所有reads根据位点的重合程度拼接在一起
     * 3. 如果有start和end site出现不止一次，那么就按照频率最高取值
     * 4. 没有特定频率取最长
     * 5. 先拼接reads本身，再拼接exon
     */
    private fun buildTemplates() {
        this.logger.info("Start to build templates")
        val geneReads = mutableMapOf<Genes, MutableList<Genes>>()

        for (i in this.matchedGeneRead) {
            if (i.gene in geneReads.keys) {
                geneReads[i.gene]!!.add(i.reads)
            } else {
                geneReads[i.gene] = mutableListOf(i.reads)
            }
        }

        val templates = mutableListOf<Template>()
        for ((k, v) in geneReads) {
            when (v.size) {
                1 -> {
                    templates.add(Template(k, v[0]))
                }
                else -> {
                    templates.add(Template(k, this.mergeReads(v)))
                }
            }
        }
        this.templates = templates
    }


    /**
     * 将两个列表的exons融合到一起
     * 就是根据上下游两个exon有没有重合位点
     * 有重合，就融合为一个
     * 没有就都保留
     * @param first 第一个列表的外显子位点
     * @param second 第二个列表的外显子位点
     * @return 全新融合有的外显子位点列表
     */
    private fun mergeExons(first: List<Array<Int>>, second: List<Array<Int>>): List<Array<Int>> {
        fun isOverlap(first: Array<Int>, second: Array<Int>): Boolean {
            return (first[0] < second[1] && first[1] > second[0])
        }

        // 生成一个列表用来收取所有的exons
        val exons = (first + second)
                .toMutableList()
                .sortedWith(compareBy( {it[0]}, {it[1]}) )

        val mergedExons = mutableListOf<Array<Int>>()

        /*
        遍历融合
        两个列表直接融合在一起排序。然后头尾遍历判断
         */
        var merged = exons[0]
        for (i in 1..(exons.size - 1)) {
            if (isOverlap(merged, exons[i])) {
                merged = arrayOf(
                        min(merged[0], exons[i][0]),
                        max(merged[1], exons[i][1])
                )
            } else {
                mergedExons.add(merged)
                merged = exons[i]
            }

            if (i == exons.size - 1) {
                mergedExons.add(merged)
            }
        }

        return mergedExons
    }


    /**
     * 用于将多条reads融合为一条的function
     * @param reads 列表，内部都是待融合的reads
     * @param freq 阈值，筛选出现频率以上的start或者end位点
     * @return 融合后的全新的reads，其实就是template的一部分了
     */
    private fun mergeReads(reads: List<Genes>, freq: Double = 60.0): Genes {
        val starts = mutableListOf<Int>()
        val ends = mutableListOf<Int>()

        for (i in reads) {
            starts.add(i.start)
            ends.add(i.end)
        }

        // 统计频率
        val start = starts.groupBy { it }
                .toList()
                .sortedByDescending { it.second.size }

        val end = ends.groupBy { it }
                .toList()
                .sortedByDescending { it.second.size }

        // 提取最终的start和end值
        val readStart = when {
            start[0].second.size / reads.size.toDouble() >= freq -> start[0].first
            else -> starts.sortedBy { it }.first()
        }

        val readEnd = when {
            end[0].second.size / reads.size.toDouble() >= freq -> end[0].first
            else -> ends.sortedBy { it }.last()
        }

        // 融合外显子
        var exons = reads.first().exons
        var i = 1
        while (i < reads.size) {
            exons = mergeExons(exons, reads[i].exons).toMutableList()
            i++
        }

        // 再手动指定一遍外显子的头和尾，保证外显子范围能够与reads范围对应
        exons.first()[0] = readStart
        exons.last()[1] = readEnd

        val newGene = Genes(
                chrom = reads.first().chrom,
                start = readStart,
                end = readEnd
        )

        newGene.exons = exons
        return newGene
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

    fun saveTemplate(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for (i in this.templates) {
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

/*
fun main(args: Array<String>) {
    val gene = GffExtractor("/home/zhang/genome/Homo_sapiens.GRCh38.91.gff3")
    val reads = BamExtractor("/home/zhang/splicehunter_test/test.bam", silent = true)
    val test = GeneReadsCoupler(gene, reads)

    test.saveTo("/home/zhang/splicehunter_test/gene_read.txt")
    test.saveTemplate("/home/zhang/splicehunter_test/templates.txt")

//    val reader = Scanner(File("/home/zhang/splicehunter_test/test.txt"))
//
//    val reads = mutableListOf<Genes>()
//    while (reader.hasNext()) {
//        val lines = reader.nextLine().split("|")[1].split("\t")
//
//        reads.add(
//                Genes(
//                        chrom=lines[0],
//                        start=lines[1].toInt(),
//                        end=lines[2].toInt()
//                )
//        )
//
//        val exons = mutableListOf<Array<Int>>()
//        val exs = lines[7].split(",")
//        for (i in 0..(exs.size - 2) step 2) {
//            exons.add(
//                    arrayOf(exs[i].toInt(), exs[i+1].toInt())
//            )
//        }
//        reads.last().exons = exons
//
//    }
//
//    reader.close()
}
*/