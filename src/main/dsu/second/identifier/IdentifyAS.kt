package dsu.second.identifier

import dsu.carrier.Exons
import dsu.carrier.GenomicLoci
import dsu.second.index.AnnotationIndex
import dsu.second.index.SJIndex
import dsu.progressbar.ProgressBar
import dsu.second.carrier.SpliceEvent
import org.apache.log4j.Logger
import kotlin.math.abs
import java.io.File
import java.io.PrintWriter

/**
 * @author Zhangyiming
 * @since 2018.09.05
 * @version 20180905
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
        private val junctions: SJIndex,
        private val reference: AnnotationIndex,
        private val overlapLevel: Double = 90.0
) {

    val logger = Logger.getLogger(IdentifyAS::class.java)
    val data: List<SpliceEvent>
    init {
        data = this.findAs(this.preProcessGjWithRef())
    }
    /**
     * 第一步将same-start和same-end的junctions
     * 按照最长的那个junctions进行配对
     *
     * 其次，先根据junctions的变化找出相应的潜在的可变剪接
     */
    fun preProcessGjWithRef(): Map<Int, List<Exons>> {
        this.logger.info("Start to paring same-start junctions with same-end junctions")

        val paired = mutableMapOf<Int, List<Exons>>()

        val pb = ProgressBar(this.junctions.sameStart.size.toLong(), "Paring same-start and same-end")

        for ( (k, v) in this.junctions.sameStart ) {
            pb.step()
            if ( v.size <= 1 ) {
                continue
            }

            val tmpValue = v.sortedByDescending { (it.end - it.start).dec() }.toMutableList()
            tmpValue.addAll(this.junctions.sameEnd[tmpValue[0].getSiteHash(false)]!!)

            // the same-start and same-end were paired by longest reads
            paired[k] = tmpValue.sorted()
        }

        return paired
    }

    private fun isSameSite(first: Int, second: Int, distanceError: Int = 0): Boolean {
        return abs(first - second) <= distanceError
    }

    /**
     * 获取最终的剪接事件
     */
    private fun getSplicingType(junctions: List<Exons>, exons: List<Exons>): List<SpliceEvent> {
        val results = mutableListOf<SpliceEvent>()

        var lastMatch = 0
        for ( i in 0..(junctions.size - 1)) {
            val currentJunction = junctions[i]
            val nextJunction = when {
                i < junctions.size - 1 -> junctions[i + 1]
                else -> currentJunction
            }
            if (
                i < junctions.size - 1 &&
                currentJunction.isOverlap(nextJunction)
            ) {
                val tmp = when {
                    currentJunction.start != nextJunction.start &&
                            currentJunction.end != nextJunction.end -> SpliceEvent("A3SS/A5SS")

                    currentJunction.start != nextJunction.start &&
                            currentJunction.end == nextJunction.end -> SpliceEvent("A5SS")

                    currentJunction.start == nextJunction.start &&
                            currentJunction.end != nextJunction.end -> SpliceEvent("A3SS")

                    else -> null
                }

                if (tmp != null) {
                    tmp.spliceIn.add(currentJunction)
                    tmp.spliceOut.add(nextJunction)
                    results.add(tmp)
                }

            }

            var j = lastMatch
            while ( j < exons.size ) {
                when {
                    currentJunction.isDownStream(exons[j], 1) -> j++
                    currentJunction.isUpStream(exons[j], 1) -> j = exons.size
                    else -> {
                        if (j != lastMatch) {
                            lastMatch = j
                        }

                        if ( currentJunction.overlapPercent(exons[j]) > this.overlapLevel ) {
                            if (
                                    currentJunction.start == exons[j].start &&
                                    currentJunction.end == exons[j].end
                            ) {
                                results.add(SpliceEvent("SE", "Exact"))
                            } else if (
                                    currentJunction.start < exons[j].start ||
                                    currentJunction.end >  exons[j].end
                                    ) {
                                results.add(SpliceEvent("SE", "Partial"))
                            } else {
                                results.add(SpliceEvent("SE", "Other"))
                            }

                            results.last().spliceIn.add(currentJunction)
                            results.last().spliceOut.add(junctions[-1])
                        }

                        if (j < exons.size - 1 && i < junctions.size) {
                            if ( currentJunction.start > exons[j].start && currentJunction.end < exons[j].end ) {
                                results.add(SpliceEvent("IE"))
                                results.last().spliceIn.add(currentJunction)
                                results.last().spliceOut.add(junctions[-1])
                            }

                            try{
                                if (
                                        GenomicLoci(
                                                currentJunction.end, nextJunction.start
                                        ).overlapPercent(
                                                GenomicLoci(exons[j].end, exons[j + 1].start)
                                        ) > this.overlapLevel
                                ) {
                                    results.add(SpliceEvent("IR"))
                                    results.last().spliceIn.addAll(listOf(currentJunction, nextJunction))
                                    results.last().spliceOut.add(junctions[-1])
                                }
                            } catch (e: dsu.errors.ChromosomeException) {

                            }

                        }

                        j ++
                    }
                }
            }
        }
        return results
    }

    /**
     * 识别所有的可变剪接
     * 将所有测序的junctions和reference中的exons排序，
     * 然后按照same-start顺序遍历一遍junctions，就可以把所有的数提取出来，
     * 以及该junctions附近的所有exons提取出来
     * @param sameStart 具有相同起始位点的Junctions
     * @param sameEnd 具有相同终止位点的Junctions
     */
    fun findAs(paired: Map<Int, List<Exons>>): List<SpliceEvent> {
        this.logger.info("Start to finding all AS")
        val tmpJunctions = this.junctions.data.keys.sorted()
        val tmpReferences = this.reference.data.sorted()
        val events = mutableListOf<SpliceEvent>()

        var i = 0
        var j = 0
        var firstMatch = 0
        while (i < tmpJunctions.size && j < tmpReferences.size) {
            val currentJunction = tmpJunctions[i]
            val currentReference = tmpReferences[j]
            when {
                currentJunction.isUpStream(currentReference) -> {
                    i++
                    j = firstMatch
                }
                currentJunction.isDownStream(currentReference) -> {
                    try{
                        val tmp = this.getSplicingType(
                                junctions = paired[currentJunction.getSiteHash(true)]!!,
                                exons = tmpReferences.subList(
                                        when (firstMatch) {
                                            0 -> 0
                                            else-> firstMatch - 1
                                        },
                                        j + 1
                                )
                        )

                        events.addAll(tmp)

                    } catch (e: NullPointerException) {

                    }
                    j++
                }

                else -> {
                    if ( j != firstMatch ) {
                        firstMatch = j
                    }

                    j ++
                }
            }
        }

        return events
    }


    fun writeTo(output: File) {
        if ( !output.absoluteFile.parentFile.exists() ) {
            output.absoluteFile.parentFile.mkdirs()
        }

        val writer = PrintWriter(output)
        for ( i in this.data ) {
            writer.println(i)
        }

        writer.close()
    }

}


fun main(args: Array<String>) {
    val pattern = "^([\\w\\.]+):(\\d+)-(\\d+)\t(\\d+)(\\+|-)$".toRegex()

    val (chromosome, start, end, strand, count) = pattern.find("ch.r1:100-200\t500-")!!.destructured

    println("$chromosome, $start, $end, $strand, $count")
}