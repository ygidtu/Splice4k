package dsu.third.template

/*
import java.io.File
import java.io.PrintWriter
import java.io.IOException
*/


import org.apache.log4j.Logger

import dsu.carrier.*
import dsu.third.extractor.*
import dsu.third.carrier.*
import dsu.progressbar.ProgressBar


/**
 * @author Zhang yiming
 * @since 2018.06.20
 * @version 20180920
 * 将基因与reads匹配到一起
 */


/**
 * 将基因与read匹配到一起的class
 * @param reference 参考基因组的Extractor
 * @param Reads 测序reads的BamExtractor
 * @param overlap 定义基因和read确实具是一对的重合程度的阈值
 * @param foldChange read与基因可能存在多对多的关系，其中重合程度最高的要大于第二多少，才将两者匹配在一起
 * @param distanceError 多少bp以内，可以认为两个位点其实是同一个位点，这里是容错率
 */
class GeneReadsCoupler(
        private val reference: Extractor,
        private val Reads: BamExtractor,
        private val overlap: Double = 90.0,
        private val foldChange: Double = 1.5,
        private val distanceError: Int = 3
        ) {
    private val logger = Logger.getLogger(GeneReadsCoupler::class.java)

    val novelReads = mutableListOf<Genes>()

    val templates = mutableListOf<Template>()

    init {
        this.matchGeneReadsWithoutChunk()
    }


    /**
     * 内存足够，不需要分块读取
     */
    private fun matchGeneReadsWithoutChunk()  {
        // 这个gap就是为了控制输出一个合适的进度条的

        val tmpMatched = mutableMapOf<Genes, Template>()
        val tmpMatchedReads = hashSetOf<Genes>()

        var firstOverlap = true
        var readIndex = 0
        this.reference.sort()
        this.Reads.sort()
        var tmpGene = this.reference.next()
        var tmpRead = this.Reads.next()

        this.logger.info("Start to matching genes and reads")
        // 统计所有的配对信息
        val pb = ProgressBar(this.reference.totalLine.toLong(), "Gene Reads matching")
        while ( tmpGene != null && tmpRead != null ) {

            pb.stepTo(this.reference.index.toLong())

            when {
                /*
                 基因在read上游，下一个基因
                 添加了3bp的误差空间，如果距离在3bp内的都算是同一个点了
                 */
                tmpGene.isUpStream(tmpRead, this.distanceError) -> {
                    tmpGene = this.reference.next()

                    // rollback read index
                    this.Reads.index = readIndex
                    tmpRead = this.Reads.get(readIndex)
                    firstOverlap = true
                }
                // 基因在read下游，读一个read
                tmpGene.isDownStream(tmpRead, this.distanceError) -> {

                    if ( tmpRead !in tmpMatchedReads ) {
                        this.novelReads.add(tmpRead)
                    }

                    tmpRead = this.Reads.next()
                }

                else -> {
                    if ( tmpGene.strand == tmpRead.strand ) {
                        // log read index
                        if (firstOverlap) {
                            readIndex = this.Reads.index
                            firstOverlap = false
                        }
                        val tmpGeneRead = GeneRead(tmpGene, tmpRead)

                        /*
                        2018.07.04
                        修正，使用基因和reads外显子的覆盖度作为基因和read匹配评判的标准

                        2018.09.05
                        这里仅用一个外显子的匹配进行配对
                        主要是为了保证能够有尽可能多的reads与reference配对，
                        在组建templates时会进行进一步的检查，保证没有错配
                         */
                        if (
                            tmpGeneRead.overlapPercent >= this.overlap &&
                            tmpGeneRead.isGeneReadsExonsOverlapQualified(this.distanceError)
                        ) {
                            // 判断是否临时的匹配中是否含有该条read了
                            if (!tmpMatched.containsKey(tmpGene)) {
                                tmpMatched[tmpGene] = Template(tmpGene, mutableListOf(tmpRead))
                            } else {
                                tmpMatched[tmpGene]!!.reads.add(tmpRead)
                            }
                            tmpMatchedReads.add(tmpRead)
                        }
                    }

                    tmpRead = this.Reads.next()
                }
            }
        }

        this.templates.addAll(tmpMatched.values)
    }

    /*
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
    */


    /*
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
    */


    /*
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
     */
}
