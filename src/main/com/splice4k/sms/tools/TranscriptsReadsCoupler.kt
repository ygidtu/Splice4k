package com.splice4k.sms.tools


import com.splice4k.base.Genes
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.sms.base.GeneRead
import com.splice4k.sms.base.Template


/**
 * @author Zhang yiming
 * @since 2018.06.20
 * @version 20180925
 * 将基因与reads匹配到一起
 *
 * 20180925 尝试将重合比例1.5倍的阈值引入，失败
 */


/**
 * 将基因与read匹配到一起的class
 * @param reference 参考基因组的Extractor
 * @param reads 测序reads的BamExtractor
 * @param overlap 定义基因和read确实具是一对的重合程度的阈值
 * @param distanceError 多少bp以内，可以认为两个位点其实是同一个位点，这里是容错率
 */
class TranscriptsReadsCoupler(
        reference: AnnotationIndex,
        reads: SJIndex,
        private val overlap: Double,
        private val distanceError: Int
        ) {

    private val novelReads = mutableListOf<Genes>()

    // 同一个基因上，转录本和reads之间的配对
    val templates = mutableMapOf<String, MutableList<Template>>()

    private val reference = reference.transcripts.sorted()
    private val reads = reads.transcripts.sorted()

    init {
        this.matchGeneReads()
        println("total: ${this.reads.size}")
        println("novel: ${this.novelReads.size}")
    }


    /**
     * 内存足够，不需要分块读取
     */
    private fun matchGeneReads()  {
        // 这个gap就是为了控制输出一个合适的进度条的

        val tmpMatched = mutableMapOf<Genes, Template>()
        val tmpMatchedReads = hashSetOf<Genes>()

        var firstOverlap = true
        var readIndex = 0

        var i = 0; var j = 0

        println("Start to matching genes and reads")
        // 统计所有的配对信息
        while ( i < this.reference.size && j < this.reads.size ) {

            val tmpGene = this.reference[i]
            val tmpRead = this.reads[j]

            when {
                /*
                 基因在read上游，下一个基因
                 添加了3bp的误差空间，如果距离在3bp内的都算是同一个点了
                 */
                tmpGene.isUpStream(tmpRead, this.distanceError) -> {
                    // rollback read index
                    if (!firstOverlap) {
                        j = readIndex
                        firstOverlap = true
                    }

                    i++
                }
                // 基因在read下游，读一个read
                tmpGene.isDownStream(tmpRead, this.distanceError) -> {

                    if ( tmpRead !in tmpMatchedReads ) {
                        this.novelReads.add(tmpRead)
                    }

                    j++
                }

                else -> {
                    if ( tmpGene.strand == tmpRead.strand ) {
                        // log read index
                        if (firstOverlap) {
                            readIndex = j
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

                    j++
                }
            }
        }

        for ( (k, v) in tmpMatched ) {
            val tmpTemplate = when ( this.templates.containsKey(k.parent) ) {
                true -> this.templates[k.parent]!!
                else -> mutableListOf()
            }

            tmpTemplate.add(v)
            this.templates[k.parent] = tmpTemplate
        }
    }
}
