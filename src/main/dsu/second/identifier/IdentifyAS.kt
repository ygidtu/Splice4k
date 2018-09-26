package dsu.second.identifier

import dsu.carrier.Exons
import dsu.carrier.GenomicLoci
import dsu.carrier.SpliceEvent
import dsu.second.index.AnnotationIndex
import dsu.second.index.SJIndex
import dsu.third.extractor.BamExtractor
import dsu.third.extractor.Extractor
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter


/**
 * @author Zhangyiming
 * @since 2018.09.05
 * @version 20180919
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


class IdentifyAS(
        val overlapOfExonIntron: Double,
        val silent: Boolean
) {
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
                        kotlin.math.abs(currentEvent.sliceSites[2] - currentExon.start) == 1 &&
                                currentEvent.sliceSites[3] in currentExon.start..currentExon.end
                ) {
                    match ++
                }
            } else if (currentEvent.sliceSites[2] == currentEvent.sliceSites[3]) {
                if (
                    kotlin.math.abs(currentEvent.sliceSites[1] - currentExon.end) == 1 &&
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

        if ( annotation == null ) {
            for ( it in events ) {
                matched[it] = mutableListOf("NA\tNA\tNA")
            }
        } else {
            var i = 0; var j = 0; var logged = 0; var firstMatch = true

            while ( i < events.size && j < annotation.size) {
                val currentEvent = events[i]
                val currentExon = annotation[j]

                when{
                    currentEvent.isUpStream( currentExon ) -> {
                        i ++

                        val qualfied = when (currentEvent.event) {
                            "SE" -> this.checkSE(currentEvent, annotation.subList(logged, j))

                            "MXE" -> this.checkMXE( currentEvent, annotation.subList(logged, j) )

                            else -> this.checkA35( currentEvent, annotation.subList(logged, j) )
                        }


                        if ( currentEvent !in matched.keys ) {
                            matched[currentEvent] = mutableListOf()
                        }

                        if ( currentEvent.sliceSites[1] ==  161036207 && currentEvent.sliceSites[2] == 161036245 && currentEvent.event == "SE") {
                            println(qualfied)
                            println(matched[currentEvent])
                        }

                        when(qualfied) {
                            null -> matched[currentEvent]!!.add("NA\tNA\tNA")
                            else -> matched[currentEvent]!!.add("${qualfied.source["gene"]}\t${qualfied.source["transcript"]}\t${qualfied.exonId}")
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
                    } catch (e: dsu.errors.ChromosomeException) {

                    }

                }
            }
        }
    }


    /**
     * 将事件与注释配对，找到事件的来源
     * @param event SJIndex
     * @param annotation Gtf或Gff注释文件index
     * @return Map of SpliceEvents and its corresponding genes, transcripts and exons
     */
    fun matchEventsWithRef( event: SJIndex, annotation: AnnotationIndex ): Map<SpliceEvent, List<String>> {
        this.logger.info("Predicting Alternative Splicing events")

        val events = mutableMapOf<String, List<SpliceEvent>>()
        val annotations = annotation.data

        for ( i in event.data.values) {
            events["${i.chromosome}${i.strand}"] = i.identifyAS(this.silent)
        }


        this.logger.info("Matching AS events with Reference")
        val res = HashMap<SpliceEvent, MutableList<String>>()

        for ( (k, v) in events ) {

            if ( k.endsWith(".") ) {
                var tmpK = k.replace("\\.$", "-")
                if ( annotations.containsKey(tmpK) ) {
                    this.matchEventsWithRefSingleChromosome(
                            v.asSequence().distinct().sorted().toList(),
                            annotations[tmpK]?.sorted(),
                            res
                    )
                }

                tmpK = k.replace("\\.$", "+")
                this.matchEventsWithRefSingleChromosome(
                        v.asSequence().distinct().sorted().toList(),
                        annotations[tmpK]?.sorted(),
                        res
                )

            } else {
                this.matchEventsWithRefSingleChromosome(
                        v.asSequence().distinct().sorted().toList(),
                        annotations[k]?.sorted(),
                        res
                )

            }
        }

        return res
    }


    /**
     * 将事件与注释配对，找到事件的来源
     * @param event SJIndex
     * @param annotation Gtf或Gff注释文件index
     * @return Map of SpliceEvents and its corresponding genes, transcripts and exons
     */
    fun matchEventsWithRef( event: BamExtractor, annotation: Extractor ): Map<SpliceEvent, List<String>> {
        this.logger.info("Predicting Alternative Splicing events")

        val events = mutableMapOf<String, List<SpliceEvent>>()
        val annotations = annotation.index

        for ( i in event.graph.values) {
            events["${i.chromosome}${i.strand}"] = i.identifyAS(this.silent)
        }


        this.logger.info("Matching AS events with Reference")
        val res = HashMap<SpliceEvent, MutableList<String>>()

        for ( (k, v) in events ) {

            if ( k.endsWith(".") ) {
                var tmpK = k.replace("\\.$", "-")
                if ( annotations.containsKey(tmpK) ) {
                    this.matchEventsWithRefSingleChromosome(
                            v.sorted(),
                            annotations[tmpK]?.sorted(),
                            res
                    )
                }

                tmpK = k.replace("\\.$", "+")
                this.matchEventsWithRefSingleChromosome(
                        v.sorted(),
                        annotations[tmpK]?.sorted(),
                        res
                )

            } else {
                this.matchEventsWithRefSingleChromosome(
                        v.sorted(),
                        annotations[k]?.sorted(),
                        res
                )

            }
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


