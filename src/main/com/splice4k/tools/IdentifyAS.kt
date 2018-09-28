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
 * @version 20180927
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


class IdentifyAS( val overlapOfExonIntron: Double ) {
    private val logger = Logger.getLogger(IdentifyAS::class.java)


    /**
     * 检查SE是否真实存在，即是否确实在注释中存在这么个外显子
     * @param currentEvent 事件
     * @param exonList 事件范围内的外显子
     * @return Exons? null -> 不匹配； Exons -> 事件相关的外显子
     */
    private fun checkSE( currentEvent: SpliceEvent, exonList: List<Exons> ): Exons? {

        for ( currentExon in exonList ) {


            if ( currentEvent.sliceSites[1] <= currentExon.start &&
                    currentEvent.sliceSites[2] >= currentExon.end ) {
                return currentExon
            }
        }
        return null
    }

    /**
     * 检查A3/A5是否真实存在，即确定是否发生在某个外显子上
     * @param currentEvent 事件
     * @param exonList 事件范围内的外显子
     * @return Exons? null -> 不匹配； Exons -> 事件相关的外显子
     */
    private fun checkA35( currentEvent: SpliceEvent, exonList: List<Exons> ): Exons? {
        var match = 0

        for ( currentExon in exonList ) {

            if ( currentEvent.sliceSites[0] == currentEvent.sliceSites[1] ) {
                if (
                        kotlin.math.abs(currentEvent.sliceSites[2] - currentExon.start) <= 1 &&
                                currentEvent.sliceSites[3] in currentExon.start..currentExon.end
                ) {
                    match ++
                }
            } else if (currentEvent.sliceSites[2] == currentEvent.sliceSites[3]) {
                if (
                    kotlin.math.abs(currentEvent.sliceSites[1] - currentExon.end) <= 1 &&
                            currentEvent.sliceSites[0] in currentExon.start..currentExon.end
                        ) {
                    match ++
                }
            }

            if ( match > 0 ) {
                return currentExon
            }
        }
        return null
    }

    /**
     * 检查MXE事件是否存在，即确定是否发生于至少两个外显子上
     * @param currentEvent 事件
     * @param exonList 事件范围内的外显子
     * @return Exons? null -> 不匹配； Exons -> 事件相关的外显子
     */
    private fun checkMXE( currentEvent: SpliceEvent, exonList: List<Exons> ): Exons? {
        var matched = 0

        for ( exon in exonList ) {
            if ( currentEvent.sliceSites[1] <= exon.start && currentEvent.sliceSites[2] >= exon.end ) {
                matched ++
            }

            if ( currentEvent.sliceSites[3] <= exon.start && currentEvent.sliceSites[4] >= exon.end ) {
                matched ++
            }

            if ( matched >= 2 ) {
                return exon
            }
        }

        return null
    }


    /**
     * 将获取到的剪接时间与基因挂钩
     */
    private fun matchEventsWithRefSingleChromosome(
            events: List<SpliceEvent>,
            annotation: List<Exons>?,
            matched: MutableMap<SpliceEvent, MutableList<String>>
    ) {

        if ( annotation != null ) {
            var i = 0; var j = 0; var logged = 0; var firstMatch = true

            while ( i < events.size && j < annotation.size) {
                val currentEvent = events[i]
                val currentExon = annotation[j]

                when{
                    currentEvent.isUpStream( currentExon ) -> {
                        i ++

                        when (currentEvent.event) {
                            "SE" -> this.checkSE(currentEvent, annotation.subList(logged, j))

                            "MXE" -> this.checkMXE( currentEvent, annotation.subList(logged, j) )

                            else -> this.checkA35( currentEvent, annotation.subList(logged, j) )
                        }?.let {
                            try {
                                if ( matched.containsKey(currentEvent) ) {
                                    matched[currentEvent]!!.add("${it.source["gene"]}\t${it.source["transcript"]}\t${it.exonId}")
                                } else {
                                    matched[currentEvent] = mutableListOf("${it.source["gene"]}\t${it.source["transcript"]}\t${it.exonId}")
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
                        val tmp1 = GenomicLoci(
                                chromosome = currentEvent.chromosome,
                                start = currentEvent.end,
                                end = events[i + 1].start
                        )

                        val tmp2 = GenomicLoci(
                                chromosome = currentExon.chromosome,
                                start = currentExon.end,
                                end = annotation[j + 1].start
                        )
                        if ( tmp2.overlapPercent(tmp1, all = true) > this.overlapOfExonIntron ) {
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

                            matched[tmp] = mutableListOf("${currentExon.source["gene"]}\t${currentExon.source["transcript"]}\t${currentExon.exonId}")
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
    ): Map<SpliceEvent, List<String>> {

        if ( show ) {
            this.logger.info("Predicting Alternative Splicing events")
        }

        val events = mutableMapOf<String, List<SpliceEvent>>()
        val res = HashMap<SpliceEvent, MutableList<String>>()

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

        if ( show ) {
            this.logger.info("Matching AS events with Reference")
        }

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


    fun writeTo(outfile: File, results: Map<SpliceEvent, List<String>>) {
        val outFile = outfile.absoluteFile

        if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

        val writer = PrintWriter(outFile)

        for ( (k, v) in results) {
            for ( j in v.distinct() ) {

                if (j != "NA\tNA\tNA") {
                    writer.println("$k\t$j")
                }
            }
        }

        writer.close()

    }
}

