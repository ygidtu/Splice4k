package dsu.carrier


import java.lang.NullPointerException
import org.apache.log4j.Logger

/**
 * @since 2018.09.19
 * @version 20190919
 * @author Zhang yiming
 *
 * 构建小型的可变剪接图，基本上就是把对应的Same Start和Same end收集到一起
 */


class SpliceGraph( val chromosome: String, val strand: Char = '+' ) {
    private val logger = Logger.getLogger(SpliceGraph::class.java)
    val starts = mutableMapOf<Int, Sites>()
    val ends = mutableMapOf<Int, Sites>()

    // 记录具有A3和A5的位点，便于最后寻找MXE
    val assSites = mapOf<String, HashSet<Int>>(
            "start" to hashSetOf(),
            "end" to hashSetOf()
    )



    /**
     * 内部类，记录所有位点信息
     */
    class Sites( val node: Int ): Comparable<Sites>, Iterator<Int> {
        val pos = mutableMapOf<Int, Int>()

        var current = 0
        var size = this.pos.size
        get() {
           return this.pos.size
        }

        /**
         * 添加同start或者同end的位点
         * @param site int 位点
         */
        fun addSite( site: Int, freq: Int? = null ) {
            if ( freq != null ) {
                this.pos[site] = freq
                return
            }

            try{
                this.pos[site] = this.pos[site]!! + 1
            } catch ( e: NullPointerException ) {
                this.pos[site] = 1
            }
        }


        /**
         * 获取离node最近的点
         */
        fun getClosestSite(): Int {
            val tmp = this.pos.keys.sorted()

            return when(
                kotlin.math.abs(tmp.first() - this.node) <
                        kotlin.math.abs(tmp.last() - this.node)
                ) {
                true -> tmp.first()
                else -> tmp.last()
            }
        }


        /**
        * 获取最长的那个位点，即距离node最远的位点
         */
        fun getExtremeSite(): Int {
            val tmp = this.pos.keys.sorted()

            return when(
                    kotlin.math.abs(tmp.first() - this.node) >
                            kotlin.math.abs(tmp.last() - this.node)
                    ) {
                true -> tmp.first()
                else -> tmp.last()
            }
        }

        /**
         * 重载
         * 便于对starts和ends进行排序
         */
        override fun compareTo(other: Sites): Int {
            return this.node - other.node
        }

        /**
         * 重载
         * 便于对Sites进行遍历
         */
        override fun next(): Int {
            this.current += 1
            return this.pos.keys.sorted()[current - 1]
        }

        override fun hasNext(): Boolean {
            return this.current < this.pos.size
        }


        fun first(): Int {
            return this.pos.keys.last()
        }

        fun last(): Int {
            return this.pos.keys.last()
        }

    }

    /**
     * 二分法进行插入
     * 保证序列有序，以及算法高效
     * @param list 包含Sites的
     */
    private fun insertSites(list: MutableMap<Int, Sites>, index: Int, site: Int, freq: Int?=null) {
        if ( list.containsKey(index) ) {
            list[index]!!.addSite(site, freq)
        } else {
            val tmp = Sites(index)
            tmp.addSite(site, freq)
            list[index] = tmp
        }
    }

    /**
     * 添加边，在图中添加junctions
     * @param start junction的start
     * @param end junction的end
     * @param startFreq start出现的次数
     * @param endFreq end出现的次数
     */
    fun addEdge(start: Int, end: Int, startFreq: Int?=null, endFreq: Int?= null) {
        insertSites(list = this.starts, index = start, site = end, freq = endFreq)
        insertSites(list = this.ends, index = end, site = start, freq = startFreq)
    }


    /**
     * 判断Same Start与Same End之间的可变剪接
     */
    private fun identifyBetweenSameStartEnd(
            starts: Sites?,
            ends: Sites?,
            res: MutableList<SpliceEvent>
    ) {

        val aS = when (this.strand) {
            '+' -> "A3"
            else -> "A5"
        }

        val aE = when (this.strand) {
            '-' -> "A5"
            else -> "A3"
        }

        val start: Int
        val end: Int

        val skippedExons = hashSetOf<Int>()
        try {
            // finding potential SE
            for ( i in starts!! ) {
                for ( j in ends!! ) {

                    if ( i < j ) {
                        res.add(SpliceEvent(
                                event = "SE",
                                chromosome = this.chromosome,
                                start = starts.node,
                                end = ends.node,
                                strand = this.strand,
                                sliceSites = listOf(starts.node, i, j, ends.node)
                        ))

                        skippedExons.add(i)
                        skippedExons.add(j)

                        this.assSites["start"]!!.add(starts.node)
                        this.assSites["end"]!!.add(ends.node)
                    }
                }
            }
        } catch (e: NullPointerException) {

        }

        // finding A3/A5
        try {
            starts!!.current = 0
            start = starts.getClosestSite()
            for ( i in starts ) {
                //
                if ( start != i && i !in skippedExons ) {
                    res.add(SpliceEvent(
                            event = aS,
                            chromosome = this.chromosome,
                            start = starts.node,
                            end = i,
                            strand = this.strand,
                            sliceSites = listOf(starts.node, starts.node, start, i)
                    ))

                    this.assSites["start"]!!.add(starts.node)
                }
            }
        } catch (e: NullPointerException) {

        }

        // finding A3/A5
        try {
            ends!!.current = 0
            end = ends.getClosestSite()
            for ( i in ends ) {
                //
                if ( end != i && i !in skippedExons ) {
                    res.add(SpliceEvent(
                            event = aE,
                            chromosome = this.chromosome,
                            start = i,
                            end = ends.node,
                            strand = this.strand,
                            sliceSites = listOf(i, end, ends.node, ends.node )
                    ))

                    this.assSites["end"]!!.add(ends.node)
                }
            }
        } catch (e: NullPointerException) {

        }
    }


    private fun isOrdered(list: List<Int>): Boolean {
        for ( i in 1..(list.size - 1)) {
            if ( list[i] < list[i - 1] ) {
                return false
            }
        }
        return true
    }

    /**
     * 找出同一条染色体上，同一条链上所有的可变剪接事件
     * @return list of SpliceEvent
     */
    fun identifyAS(): List<SpliceEvent> {

        this.logger.info("Finding Alternative splicing events of ${this.chromosome}${this.strand}")
        val res = mutableListOf<SpliceEvent>()

        val loggedEnds = HashSet<Sites>()

        // 识别可变剪接
        for ( i in this.starts.values.sorted() ) {
            if (i.size <= 1) {
                continue
            }

            val tmpEnds = this.ends[i.getExtremeSite()]
            identifyBetweenSameStartEnd(
                    starts = i,
                    ends = tmpEnds,
                    res = res
            )

            if ( tmpEnds != null ) {
                loggedEnds.add(tmpEnds)
            }
        }

        for ( i in this.ends.values.sorted() ) {
            if (i.size <= 1) {
                continue
            }

            if ( !loggedEnds.contains(i) ) {
                identifyBetweenSameStartEnd(
                        starts = null,
                        ends = i,
                        res = res
                )
            }
        }

        // finding MXE
        var sites = this.starts.keys.asSequence().sorted().toMutableList()
        sites.addAll(this.ends.keys)
        sites = sites.asSequence().distinct().sorted().toMutableList()

        val assStarts = this.assSites["start"]!!.asSequence().sorted().distinct().toList()
        val assEnds = this.assSites["end"]!!.asSequence().distinct().sorted().toList()

        var firstEnds = 0
        for ( i in 0..(assStarts.size - 1)) {
            var j = firstEnds

            while ( j < assEnds.size - 1 ) {

                // 如果j太小，就跳过
                if ( assEnds[j] <= this.starts[assStarts[i]]!!.getExtremeSite() ) {
                    if ( firstEnds != j ) {

                    }
                    j ++
                    firstEnds = j
                    continue
                }

                val currentStart = this.starts[assStarts[i]]!!.pos.keys.toList()
                val currentEnd = this.ends[assEnds[j]]!!.pos.keys.toList()

                for ( e1 in 0..(currentStart.size - 1) ) {
                    for ( e2 in (e1 + 1)..(currentStart.size - 1) ) {
                        for ( s1 in 0..(currentEnd.size - 1) ) {
                            for ( s2 in (s1 + 1)..(currentEnd.size - 1) ) {
                                val mxeSites = listOf(
                                        assStarts[i],
                                        currentStart[e1],
                                        currentEnd[s1],
                                        currentStart[e2],
                                        currentEnd[s2],
                                        assEnds[j]
                                )

                                if (
                                        this.isOrdered(mxeSites) &&
                                        sites.indexOf(currentStart[e2]) - sites.indexOf(currentEnd[s1]) == 1
                                ) {
                                    res.add(SpliceEvent(
                                            event = "MXE",
                                            chromosome = this.chromosome,
                                            start = assStarts[i],
                                            end = assEnds[j],
                                            strand = this.strand,
                                            sliceSites = mxeSites
                                    ))
                                }

                            }
                        }

                    }

                }

                j++

            }
        }

        return res
    }
}