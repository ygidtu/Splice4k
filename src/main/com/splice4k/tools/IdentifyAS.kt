package com.splice4k.tools

import com.splice4k.base.Exons
import com.splice4k.base.GenomicLoci
import com.splice4k.base.JunctionsGraph
import com.splice4k.base.SpliceEvent
import com.splice4k.errors.ChromosomeException
import com.splice4k.index.AnnotationIndex
import org.apache.log4j.Logger
import java.io.File
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.Future


/**
 * @author Zhangyiming
 * @since 2018.09.05
 * @version 20180929
 *
 * 这里主要是根据二代测序的可变剪接识别方式进行AS的识别
 *
 * 第一步将same start和same end的junctions
 * 按照最长的那个junctions进行配对
 * 配对完成后可以先根据junctions的变化找出相应的可变剪接（eg: A3SS, A5SS, potential SE and IR）
 *
 * 第二步将reference中的exon与junctions进行对比，确定每个剪接的来源转录本，和对最终的可变剪接进行定性
 *
 * 第三步计算psi
 *
 */


/**
 * @param overlapOfExonIntron 识别IR所需的外显子与内含子的重合程度
 * @param bamFile 提供用来找IR PSI的bam文件
 */
class IdentifyAS(
        private val overlapOfExonIntron: Double,
        private val bamFile: File?
) {
    private val logger = Logger.getLogger(IdentifyAS::class.java)
    private val check = CheckAS()


    /**
     * 为多进程写的类，继承自Callable
     * @param graph 构建的same start和same end图
     * @param bamFile 提供用来找IR PSI的bam文件
     * @param annotations 参考基因组的基因（转录本）以及对应的外显子
     * @param overlapOfExonIntron 判断IR时外显子与内含子之间重合程度
     * @param error 判断A3/A5是否存（等）时所需要的阈值
     * @param show 是否输出详细信息
     * @param logger 日志
     */
    class Run(
            private val graph: JunctionsGraph,
            private val bamFile: File?,
            private val annotations: Map<String, List<Exons>>,
            private val overlapOfExonIntron: Double,
            private val error: Int,
            private val show: Boolean,
            private val logger: Logger
    ): Callable<MutableMap<SpliceEvent, MutableList<Exons>>> {
        private val psiOfIR = PsiOfIR()
        private val checkAS = CheckAS()

        /**
         * 将获取到的剪接时间与基因挂钩
         */
        private fun matchEventsWithRefSingleChromosome(
                events: List<SpliceEvent>,
                annotation: List<Exons>?
        ): MutableMap<SpliceEvent, MutableList<Exons>> {

            // 事件以及匹配到一起的外显子
            val matched = mutableMapOf<SpliceEvent, MutableList<Exons>>()

            if ( annotation != null ) {
                var i = 0; var j = 0; var logged = 0; var firstMatch = true

                while ( i < events.size && j < annotation.size) {
                    val currentEvent = events[i]
                    val currentExon = annotation[j]

                    when{
                        currentEvent.isUpStream( currentExon ) -> {
                            i ++

                            this.checkAS.check( currentEvent, annotation.subList(logged, j) )?.let {
                                try {

                                    if ( matched.containsKey(currentEvent) ) {
                                        matched[currentEvent]!!.addAll( it )
                                    } else {
                                        matched[currentEvent] = it.toMutableList()
                                    }

                                } catch (error: kotlin.KotlinNullPointerException ) {
                                    this.logger.error(error)

                                    for ( e in error.stackTrace ) {
                                        this.logger.error(e)
                                    }
                                    this.logger.error("Event is $currentEvent")
                                }
                            }

                            if ( !firstMatch ) {
                                j = logged
                                firstMatch = true
                            }
                        }
                        currentEvent.isDownStream( currentExon ) -> {
                            j++
                        }
                        else -> {

                            if ( firstMatch ) {
                                logged = j
                                firstMatch = false
                            }

                            j++
                        }
                    }

                    if (
                            i < events.size - 1 &&
                            j < annotation.size - 1 &&
                            currentExon.source == annotation[j + 1].source
                    ) {
                        try{
                            // junction gap
                            val tmp1 = GenomicLoci(
                                    chromosome = currentEvent.chromosome,
                                    start = currentEvent.end,
                                    end = events[i + 1].start
                            )

                            // exons gap
                            val tmp2 = GenomicLoci(
                                    chromosome = currentExon.chromosome,
                                    start = currentExon.end,
                                    end = annotation[j + 1].start
                            )
                            if (

                                    tmp2.overlapPercent(tmp1, all = true) > this.overlapOfExonIntron
                            ) {
                                val tmp = SpliceEvent(
                                        event = "IR",
                                        chromosome = currentEvent.chromosome,
                                        start = kotlin.math.min(currentEvent.end, currentExon.end),
                                        end = kotlin.math.max(events[i + 1].start, annotation[j + 1].start),
                                        strand = currentEvent.strand,
                                        sliceSites = mutableListOf(
                                                currentEvent.end,
                                                currentExon.end,
                                                events[i + 1].start,
                                                annotation[j + 1].start
                                        )
                                )
                                tmp.isNovel = false
                                this.bamFile?.let {
                                    tmp.psi = this.psiOfIR.getPsi(
                                            chromosome = currentEvent.chromosome,
                                            regionStart = currentExon.end,
                                            regionEnd = annotation[j + 1].start,
                                            bamFile = this.bamFile
                                    )
                                }
                                matched[tmp] = mutableListOf(currentExon)
                            }
                        } catch (e: ChromosomeException) {

                        }

                    }
                }
            }

            return matched
        }


        /**
         * 重载call function
         * @return 与外部class的结果类别一致
         */
        override fun call(): MutableMap<SpliceEvent, MutableList<Exons>>  {
            val tmpRes = this.graph.identifyAS(error = error, silent = !show)
            val k = "${this.graph.chromosome}${this.graph.strand}"

            val results = mutableMapOf<SpliceEvent, MutableList<Exons>>()
            if ( k.endsWith(".") ) {
                var tmpK = k.replace("\\.$", "-")
                if ( annotations.containsKey(tmpK) ) {
                    results.putAll(
                            this.matchEventsWithRefSingleChromosome(
                                    events = tmpRes.asSequence().distinct().sorted().toList(),
                                    annotation = annotations[tmpK]?.sorted()
                            )
                    )
                }

                tmpK = k.replace("\\.$", "+")
                results.putAll(
                        this.matchEventsWithRefSingleChromosome(
                                events = tmpRes.asSequence().distinct().sorted().toList(),
                                annotation = annotations[tmpK]?.sorted()
                        )
                )
            } else {
                results.putAll(
                        this.matchEventsWithRefSingleChromosome(
                                events = tmpRes.asSequence().distinct().sorted().toList(),
                                annotation = annotations[k]?.sorted()
                        )
                )
            }
            return results
        }
    }


    /**
     * 将事件与注释配对，找到事件的来源
     * @param event List of SpliceGraph
     * @param annotations Gtf或Gff注释文件提取的exon index， Map of gene, transcript and its corresponding exons
     * @return Map of SpliceEvents and its corresponding genes, transcripts and exons
     */
    fun matchEventsWithRef(
            event: List<JunctionsGraph>,
            annotations: AnnotationIndex,
            threads: Int,
            error: Int,
            show: Boolean = true
    ): Map<SpliceEvent, List<Exons>> {

        this.logger.info("Predicting Alternative Splicing events")

        // val events = mutableMapOf<String, List<SpliceEvent>>()
        val res = HashMap<SpliceEvent, MutableList<Exons>>()
        val pool = Executors.newFixedThreadPool( threads )
        val futures = mutableListOf<Future<MutableMap<SpliceEvent, MutableList<Exons>>>>()


        for ( i in event ) {
            if ( i.isEmpty() ) {
                continue
            }

            val f = pool.submit(Run(
                    graph = i,
                    bamFile = this.bamFile,
                    annotations = annotations.data,
                    overlapOfExonIntron = this.overlapOfExonIntron,
                    error = error,
                    show = show,
                    logger = this.logger
            ))

            futures.add(f)
        }

        futures.forEach {
            res.putAll(it.get())
        }

        pool.shutdown()

        check.checkALEAFE(res, annotations.transcripts)

        return res
    }
}


