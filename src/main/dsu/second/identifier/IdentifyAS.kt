package dsu.second.identifier


import dsu.carrier.Exons
import dsu.carrier.GenomicLoci
import dsu.second.index.AnnotationIndex
import dsu.second.index.SJIndex
import dsu.progressbar.ProgressBar
import dsu.carrier.SpliceEvent
import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
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


class IdentifyAS(event: SJIndex, annotation: AnnotationIndex) {
    private val logger = Logger.getLogger(IdentifyAS::class.java)
    private val events = mutableMapOf<String, List<SpliceEvent>>()
    private val annotation = annotation.data

    init {
        for ( i in event.data.values) {
            this.events["${i.chromosome}${i.strand}"] = i.identifyAS()
        }
    }

    /**
     * 将获取到的剪接时间与基因挂钩
     */
    fun matchEventsWithRefSingleChromosome(
            events: List<SpliceEvent>,
            annotation: List<Exons>
    ): Map<SpliceEvent, String> {
        var i = 0
        var j = 0

        val matched = mutableMapOf<SpliceEvent, String>()
        while ( i < events.size && j < annotation.size) {
            val currentEvent = events[i]
            val currentExon = annotation[j]

            when{
                currentEvent.isUpStream( currentExon, distanceError = 3 ) -> {
                    if ( !matched.containsKey(currentEvent) ) {
                        matched[currentEvent] = ""
                    }
                    i ++
                }
                currentEvent.isDownStream( currentEvent, distanceError = 3 ) -> j++
                else -> {
                    if ( currentEvent.strand == currentExon.strand ) {
                        matched[currentEvent] = currentExon.source
                    }
                    j++
                }

            }

            if (
                    i < this.events.size - 1 &&
                    j < this.annotation.size &&
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
                    if ( tmp1.overlapPercent(tmp2) > 90.0 ) {
                        val tmp = SpliceEvent(
                                event = "IR",
                                chromosome = currentEvent.chromosome,
                                start = kotlin.math.min(currentEvent.end, currentExon.end),
                                end = kotlin.math.max(events[i + 1].start, annotation[j + 1].start),
                                strand = currentEvent.strand,
                                sliceSites = listOf(
                                        currentEvent.end,
                                        currentExon.end,
                                        events[i + 1].start,
                                        annotation[j + 1].start
                                )
                        )

                        matched[tmp] = currentExon.source
                    }
                } catch (e: dsu.errors.ChromosomeException) {

                }

            }
        }

        return matched
    }


    fun matchEventsWithRef(): Map<SpliceEvent, String> {
        val res = HashMap<SpliceEvent, String>()
        for ( (k, v) in this.events ) {
            if ( this.annotation.containsKey(k) ) {
                res.putAll(this.matchEventsWithRefSingleChromosome(v, this.annotation[k]!!))
            }
        }

        return res
    }


    fun writeTo(outfile: File) {
        val outFile = outfile.absoluteFile

        var writer = PrintWriter(System.out)

        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for ( (k, v) in this.matchEventsWithRef()) {
                writer.println("$k\t$v")
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


