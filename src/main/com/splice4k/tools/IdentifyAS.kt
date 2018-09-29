package com.splice4k.tools

import com.splice4k.base.Exons
import com.splice4k.base.GenomicLoci
import com.splice4k.base.SpliceEvent
import com.splice4k.base.SpliceGraph
import com.splice4k.errors.ChromosomeException
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter
import java.util.concurrent.Executors



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
 */
class IdentifyAS(
        private val overlapOfExonIntron: Double,
        private val bamFile: File?
) {
    private val logger = Logger.getLogger(IdentifyAS::class.java)
    private val psiOfIR = PsiOfIR()
    private val checkAS = CheckAS()


    /**
     * 将获取到的剪接时间与基因挂钩
     */
    private fun matchEventsWithRefSingleChromosome(
            events: List<SpliceEvent>,
            annotation: List<Exons>?,
            matched: MutableMap<SpliceEvent, MutableList<Exons>>
    ) {

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
                                    matched[currentEvent]!!.add( it )
                                } else {
                                    matched[currentEvent] = mutableListOf( it )
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
    }


    /**
     * 将事件与注释配对，找到事件的来源
     * @param event List of SpliceGraph
     * @param annotations Gtf或Gff注释文件提取的exon index， Map of gene, transcript and its corresponding exons
     * @return Map of SpliceEvents and its corresponding genes, transcripts and exons
     */
    fun matchEventsWithRef(
            event: List<SpliceGraph>,
            annotations: Map<String, List<Exons>>,
            threads: Int,
            error: Int,
            show: Boolean = true
    ): Map<SpliceEvent, List<Exons>> {

        this.logger.info("Predicting Alternative Splicing events")

        val events = mutableMapOf<String, List<SpliceEvent>>()
        val res = HashMap<SpliceEvent, MutableList<Exons>>()

        var executor = Executors.newFixedThreadPool(threads)

        for ( i in event) {
            val worker = Runnable {
                val tmpEvents = i.identifyAS(error = error, silent = !show)
                events["${i.chromosome}${i.strand}"] = tmpEvents
            }
            executor.execute(worker)
        }
        executor.shutdown()
        while (!executor.isTerminated) {

        }

        executor = Executors.newFixedThreadPool(threads)

        this.logger.info("Matching AS events with Reference")

        for ( (k, v) in events ) {

            if ( k.endsWith(".") ) {
                var tmpK = k.replace("\\.$", "-")
                if ( annotations.containsKey(tmpK) ) {

                    val worker = Runnable {
                        this.matchEventsWithRefSingleChromosome(
                                events = v.asSequence().distinct().sorted().toList(),
                                annotation = annotations[tmpK]?.sorted(),
                                matched = res
                        )
                    }
                    executor.execute(worker)
                }

                tmpK = k.replace("\\.$", "+")
                val worker = Runnable {
                    this.matchEventsWithRefSingleChromosome(
                            events = v.asSequence().distinct().sorted().toList(),
                            annotation = annotations[tmpK]?.sorted(),
                            matched = res
                    )
                }
                executor.execute(worker)

            } else {
                val worker = Runnable {
                    this.matchEventsWithRefSingleChromosome(
                            events = v.asSequence().distinct().sorted().toList(),
                            annotation = annotations[k]?.sorted(),
                            matched = res
                    )
                }
                executor.execute(worker)

            }
        }


        executor.shutdown()
        while (!executor.isTerminated) {
        }

        return res
    }

    fun writeTo(outfile: File, results: Map<SpliceEvent, List<Exons>>) {
        val outFile = outfile.absoluteFile

        if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

        val writer = PrintWriter(outFile)

        writer.println("#spliceRange\tspliceType\tspliceSites\tgene\ttranscript\texon\tPSI")
        for ( (k, v) in results) {
            for ( j in v.distinct() ) {
                writer.println("$k\t$j\t${k.psi}")
            }
        }

        writer.close()

    }
}


