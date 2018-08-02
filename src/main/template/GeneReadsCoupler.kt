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


/**
 * 将基因与read匹配到一起的class
 * @param Gene 参考基因组的Extractor
 * @param Reads 测序reads的BamExtractor
 * @param overlap 定义基因和read确实具是一对的重合程度的阈值
 * @param foldChange read与基因可能存在多对多的关系，其中重合程度最高的要大于第二多少，才将两者匹配在一起
 * @param distanceError 多少bp以内，可以认为两个位点其实是同一个位点，这里是容错率
 * @param chunk 若内存不够，文件过大，分成小块读取
 */
class GeneReadsCoupler(
        private val Gene: Extractor,
        private val Reads: BamExtractor,
        private val overlap: Double = 90.0,
        private val foldChange: Double = 1.5,
        private val distanceError: Int = 3,
        private val chunk: Boolean = false
        ) {
    private val logger = Logger.getLogger(GeneReadsCoupler::class.java)

    val matchedGeneRead = mutableListOf<GeneRead>()
    val novelReads = mutableListOf<Genes>()
    val fusionReads = mutableListOf<GeneRead>()

    lateinit var templates : MutableList<Template>

    init {
        this.matchGeneReads()
    }


    /**
     * 内存足够，不需要分块读取
     */
    private fun matchGeneReadsWithoutChunk(): MutableMap<Genes, MutableList<GeneRead>>  {
        // 这个gap就是为了控制输出一个合适的进度条的
        val gap = this.Gene.totalLine.toString().length - 2
        var logged = -1
        val tmpMatched = mutableMapOf<Genes, MutableList<GeneRead>>()

        var firstOverlap = true
        var readIndex = 0
        var tmpGene = this.Gene.next()
        var tmpRead = this.Reads.next()

        // 统计所有的配对信息
        while (tmpGene != null && tmpRead != null) {
            if (this.Gene.index % 10.0.pow(gap.toDouble()).toInt() == 0) {
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
                    if (!tmpMatched.containsKey(tmpRead)) {
                        this.novelReads.add(tmpRead)
                    }

                    tmpRead = this.Reads.next()
                }

                else -> {
                    if (tmpGene.strand == tmpRead.strand) {
                        // log read index
                        if (firstOverlap) {
                            readIndex = this.Reads.index
                            firstOverlap = false
                        }
                        val tmpGeneRead = GeneRead(tmpGene, tmpRead)

                        /*
                        2018.07.04
                        修正，使用基因和reads外显子的覆盖度作为基因和read匹配评判的标准
                         */
                        if (tmpGeneRead.isGeneReadsExonsOverlapQualified(this.overlap)) {
                            // 判断是否临时的匹配中是否含有该条read了
                            if (tmpMatched.containsKey(tmpRead)) {
                                tmpMatched[tmpRead]!!.add(tmpGeneRead)
                            } else {
                                tmpMatched[tmpRead] = mutableListOf(tmpGeneRead)
                            }
                        }
                    }

                    tmpRead = this.Reads.next()
                }
            }
        }
        return tmpMatched
    }


    /**
     * 文件过大，内存不够，分割文件进行读取
     * 唯独一点需要注意，为了保证有3bp的容错，
     * 在获取overlap区域的时候，需要在基因的start和end区域进行相应的加减
     * 由于这个获取方式，可能会导致不停的遍历整个bam文件
     * 因此执行效率可能会有比较大的降低
     */
    private fun matchGeneReadsWithChunk(): MutableMap<Genes, MutableList<GeneRead>> {
        val tmpMatched = mutableMapOf<Genes, MutableList<GeneRead>>()

        val gap = this.Gene.totalLine.toString().length - 2
        var logged = -1

        var tmpGene = this.Gene.next()

        while (tmpGene != null) {
            if (this.Gene.index % 10.0.pow(gap.toDouble()).toInt() == 0) {
                if (logged != this.Gene.index) {
                    logger.info("Gene Reads matching at ${this.Gene.index}/${this.Gene.totalLine}")
                    logged = this.Gene.index
                }
            }

            val tmpReads = this.Reads.getAllRegions(
                    tmpGene.chrom,
                    tmpGene.start - this.distanceError,
                    tmpGene.end + this.distanceError
            )


            for (tmpRead in tmpReads) {
                if (tmpGene.strand == tmpRead.strand) {
                    val tmpGeneRead = GeneRead(tmpGene, tmpRead)

                    /*
                    2018.07.04
                    修正，使用基因和reads外显子的覆盖度作为基因和read匹配评判的标准
                     */
                    if (tmpGeneRead.isGeneReadsExonsOverlapQualified(this.overlap)) {
                        // 判断是否临时的匹配中是否含有该条read了
                        if (tmpMatched.containsKey(tmpRead)) {
                            tmpMatched[tmpRead]!!.add(tmpGeneRead)
                        } else {
                            tmpMatched[tmpRead] = mutableListOf(tmpGeneRead)
                        }
                    }
                }
            }

            tmpGene = this.Gene.next()
        }

        return tmpMatched
    }


    /**
     * 将基因与reads匹配到一起
     * 根据是否采用chunk，
     * 采用不同的文件读取方式进行处理
     */
    private fun matchGeneReads() {
        val tmpMatched = when {
            this.chunk -> this.matchGeneReadsWithoutChunk()
            else -> this.matchGeneReadsWithChunk()
        }

        // 添加重合程度判断，相较于两者中短的区域，覆盖程度要大于90%
        for (v in tmpMatched.values) {
            when {
                v.size > 1 -> {
                    val tmpV = v.sortedBy { it.overlap.dec() }
                    val overlapFC = tmpV[0].overlap / tmpV[1].overlap.toDouble()

                    if ( overlapFC > this.foldChange ) {
                        if (tmpV[0].overlapPercent > this.overlap) this.matchedGeneRead.add(tmpV[0])
                    } else {
                        this.fusionReads.addAll(tmpV)
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

        this.fusionReads.sortWith(compareBy({it.reads}, {it.gene}))

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

    /**
     * 保存novel的read到文件
     * @param outfile 输出文件路径
     */
    fun saveNovel(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

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

    /**
     * 保存构件好的基因的template
     * @param outfile 输出文件路径
     */
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


    /**
     * 保存可能的融合基因
     * @param outfile 输出文件路径
     */
    fun savePotentialFusion(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for (i in this.fusionReads) {
                writer.println(i.toStringReverse())
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
