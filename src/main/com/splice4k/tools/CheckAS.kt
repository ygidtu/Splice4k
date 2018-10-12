package com.splice4k.tools

import com.splice4k.base.Exons
import com.splice4k.base.SpliceEvent

/**
 * 检查各种剪接时间发生于哪个外显子上
 */

class CheckAS() {
    /**
     * 检查SE是否真实存在，即是否确实在注释中存在这么个外显子
     * @param currentEvent 事件
     * @param exonList 事件范围内的外显子
     * @return Exons? null -> 不匹配； Exons -> 事件相关的外显子
     */
    private fun checkSE(currentEvent: SpliceEvent, exonList: List<Exons> ): List<Exons> {
        val res = mutableListOf<Exons>()
        for ( currentExon in exonList ) {

            if ( currentEvent.sliceSites[1] <= currentExon.start &&
                    currentEvent.sliceSites[2] >= currentExon.end ) {

                currentEvent.isNovel = false
                res.add(currentExon)
            }
        }
        return res
    }

    /**
     * 检查A3/A5是否真实存在，即确定是否发生在某个外显子上
     * @param currentEvent 事件
     * @param exonList 事件范围内的外显子
     * @return Exons? null -> 不匹配； Exons -> 事件相关的外显子
     */
    private fun checkA35(currentEvent: SpliceEvent, exonList: List<Exons> ): List<Exons> {
        val match = mutableSetOf<Exons>()

        for ( currentExon in exonList ) {

            if ( currentEvent.sliceSites[0] == currentEvent.sliceSites[1] ) {
                if (
                        kotlin.math.abs(currentEvent.sliceSites[2] - currentExon.start) <= 3 ||
                        kotlin.math.abs(currentEvent.sliceSites[3] - currentExon.start) <= 3
                ) {
                    match.add(currentExon)
                }
            } else if (currentEvent.sliceSites[2] == currentEvent.sliceSites[3]) {
                if (
                        kotlin.math.abs(currentEvent.sliceSites[1] - currentExon.end) <= 3 ||
                        kotlin.math.abs(currentEvent.sliceSites[1] - currentExon.end) <= 3
                ) {
                    match.add(currentExon)
                }
            }
        }

        if ( match.size > 1 ) {
            currentEvent.isNovel = false
        }

        return match.toList()
    }

    /**
     * 检查MXE事件是否存在，即确定是否发生于至少两个外显子上
     * @param currentEvent 事件
     * @param exonList 事件范围内的外显子
     * @return Exons? null -> 不匹配； Exons -> 事件相关的外显子
     */
    private fun checkMXE(currentEvent: SpliceEvent, exonList: List<Exons> ): List<Exons> {
        val matched1 = mutableSetOf<Exons>()
        val matched2 = mutableSetOf<Exons>()
        val geneExons = mutableMapOf<String, MutableSet<Exons>>()

        for ( exon in exonList ) {
            if ( currentEvent.sliceSites[1] <= exon.start && currentEvent.sliceSites[2] >= exon.end ) {
                matched1.add(exon)

                exon.source["gene"]!!.forEach {
                    val tmp = mutableSetOf(exon)
                    if ( geneExons.containsKey(it) ) {
                        tmp.addAll(geneExons[it]!!)
                    }
                    geneExons[it] = tmp
                }
            }

            if ( currentEvent.sliceSites[3] <= exon.start && currentEvent.sliceSites[4] >= exon.end ) {
                matched2.add(exon)

                exon.source["gene"]!!.forEach {
                    val tmp = mutableSetOf(exon)
                    if ( geneExons.containsKey(it) ) {
                        tmp.addAll(geneExons[it]!!)
                    }
                    geneExons[it] = tmp
                }
            }
        }


        /*
         稍微写一下思路，有点复杂
         已经在上边的遍历收集了所有的基因及其外现在
         那么我们就来看，如果同一个基因有两个及以上的外显子，
         且这两个外显子分别能够跟MXE的两个skipped区域对应，
         那么就认为这个MXE来自这个基因
          */
        if ( matched1.isNotEmpty() && matched2.isNotEmpty() ) {

            for ( ( k, values ) in geneExons ) {
                if (values.size >= 2) {
                    var match = 0
                    val matched = mutableListOf<Exons>()
                    for ( v in values.sorted()) {
                        if ( v in matched1 ) {
                            match ++
                            matched.add(v)
                        }

                        if ( v in matched2 ) {
                            match *= -1
                            matched.add(v)
                        }

                        if ( match < 0 ) {
                            val res = Exons(
                                    chromosome = v.chromosome,
                                    start = v.start,
                                    end = v.end,
                                    exonId = matched.asSequence().map { it.exonId }.joinToString(prefix = "", postfix = "", separator = ",")
                            )

                            res.source["gene"]!!.add(k)

                            matched.forEach { res.source["transcript"]!!.addAll(it.source["transcript"]!!) }

                            currentEvent.isNovel = false
                            return listOf(res)
                        }
                    }
                }
            }
        }

        return when (matched1.isNotEmpty()) {
            true -> matched1.toMutableList()
            else -> matched2.toMutableList()
        }
    }


    fun check( currentEvent: SpliceEvent, exonList: List<Exons> ): List<Exons>? {
        return when ( currentEvent.event ) {
            "A3" -> this.checkA35(currentEvent, exonList)
            "A5" -> this.checkA35(currentEvent, exonList)
            "SE" -> this.checkSE(currentEvent, exonList)
            "MXE" -> this.checkMXE(currentEvent, exonList)
            else -> null
        }
    }
}