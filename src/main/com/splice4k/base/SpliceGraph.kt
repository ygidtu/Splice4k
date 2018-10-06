package com.splice4k.base


import org.apache.log4j.Logger
import java.util.Objects
import kotlin.collections.HashSet
import kotlin.collections.component1
import kotlin.collections.component2
import kotlin.collections.set
import kotlin.math.abs


/**
 * @since 2018.09.19
 * @version 20191006
 * @author Zhang yiming
 *
 * 构建小型的可变剪接图，基本上就是把对应的Same Start和Same end收集到一起
 */


/**
 * @param chromosome 染色体
 * @param strand 链
 */
class SpliceGraph(
        val chromosome: String,
        val strand: Char
) {

    private val logger = Logger.getLogger(SpliceGraph::class.java)
    private val starts = mutableMapOf<Int, Sites>()
    private val ends = mutableMapOf<Int, Sites>()

    // 记录具有A3和A5的位点，便于最后寻找MXE
    private val as35Sites = mapOf<String, HashSet<Int>>(
            "start" to hashSetOf(),
            "end" to hashSetOf()
    )


    /**
     * 二分法进行插入
     * 保证序列有序，以及算法高效
     * @param list 包含Sites的list，主要为class内的starts和ends
     * @param index 指每个same starts和same ends的node
     * @param site 指每个same starts上对应的ends，与每个same ends上对应的starts
     * @param freq 出现的频率，由于总是成对导入，因此一个freq就够了
     * @param transcript 该starts或ends所在的转录本
     * @param gene 该starts或ends所在的基因
     */
    private fun insertSites(
            list: MutableMap<Int, Sites>,
            index: Int,
            site: Int,
            freq: Int?=null,
            transcript: String? = null,
            gene: String? = null
    ) {
        if ( list.containsKey(index) ) {
            list[index]!!.addSite(
                    site = site,
                    freq = freq,
                    transcript = transcript,
                    gene = gene
            )
        } else {
            val tmp = Sites(index)
            tmp.addSite(
                    site = site,
                    freq = freq,
                    transcript = transcript,
                    gene = gene
            )
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
    fun addEdge(
            start: Int,
            end: Int,
            freq: Int?=null,
            transcript: String? = null,
            gene: String? = null
    ) {
        this.insertSites(
                list = this.starts,
                index = start,
                site = end,
                freq = freq,
                transcript = transcript,
                gene = gene
        )
        this.insertSites(
                list = this.ends,
                index = end,
                site = start,
                freq = freq,
                transcript = transcript,
                gene = gene
        )
    }


    /**
     * 判断Same Start与Same End之间的可变剪接
     */
    private fun identifyBetweenSameStartEnd(
            starts: Sites?,
            ends: Sites?,
            res: MutableList<SpliceEvent>,
            error: Int
    ) {
        val exonSkipped = hashSetOf<Int>()

        try {
            // finding potential SE
            for ( i in starts!!.getSites() ) {
                for ( j in ends!!.getSites() ) {

                    if ( i < j ) {

                        val tmpEvent = SpliceEvent(
                                event = "SE",
                                chromosome = this.chromosome,
                                start = i.site,
                                end = j.site,
                                strand = this.strand,
                                sliceSites = mutableListOf(starts.node, i.site, j.site, ends.node)
                        )

                        tmpEvent.psi = (starts.getPsi(i.site)  + ends.getPsi(j.site)) / 2

                        res.add(tmpEvent)

                        exonSkipped.add(Objects.hash(listOf(starts.node, i.site, j.site).sorted()))
                        exonSkipped.add(Objects.hash(listOf(i.site, j.site, ends.node)))

                        this.as35Sites["start"]!!.add(starts.node)
                        this.as35Sites["end"]!!.add(ends.node)
                    }
                }
            }
        } catch (e: NullPointerException) {

        }

        try {
            // finding A3/A5
            for ( i in 0..(starts!!.size - 2) ) {
                for ( j in (i + 1)..(starts.size - 1) ) {
                    val sites = mutableListOf(
                            starts.node,
                            starts.node,
                            starts.getSite(i)!!.site,
                            starts.getSite(j)!!.site
                    )

                    sites.sort()
                    if (
                            abs(starts.getSite(i)!!.site - starts.getSite(j)!!.site) >= error &&
                            Objects.hash(sites.asSequence().sorted().distinct()) !in exonSkipped
                    ) {

                        val tmpEvent = SpliceEvent(
                                event = when(this.strand) {
                                    '-' -> "A5"
                                    else -> "A3"
                                },
                                chromosome = this.chromosome,
                                start = sites[2],
                                end = sites[3],
                                strand = this.strand,
                                sliceSites = sites
                        )

                        tmpEvent.psi = starts.getPsi( target = sites[2] )
                        res.add(tmpEvent)

                        this.as35Sites["start"]!!.add(starts.node)
                    }
                }
            }
        } catch (e: NullPointerException) {

        }

        try{
            // finding A3/A5
            for ( i in 0..(ends!!.size - 2) ) {
                for ( j in (i + 1)..(ends.size - 1) ) {

                    val sites = mutableListOf(
                            ends.getSite(i)!!.site,
                            ends.getSite(j)!!.site,
                            ends.node,
                            ends.node
                    )
                    sites.sort()
                    if (
                            abs(ends.getSite(i)!!.site - ends.getSite(j)!!.site) > error &&
                            Objects.hash(sites.asSequence().sorted().distinct()) !in exonSkipped
                    ) {

                        val tmpEvent = SpliceEvent(
                                event = when(this.strand) {
                                    '-' -> "A3"
                                    else -> "A5"
                                },
                                chromosome = this.chromosome,
                                start = sites[0],
                                end = sites[1],
                                strand = this.strand,
                                sliceSites = sites
                        )

                        tmpEvent.psi = ends.getPsi( target = sites[1] )

                        res.add(tmpEvent)

                        this.as35Sites["end"]!!.add(ends.node)
                    }
                }
            }
        } catch (e: NullPointerException) {

        }

    }


    /**
     * 判断一个列表是否已经排好序
     * @param list input list
     * @return bool true -> sorted; false -> not
     */
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
     * @param error 识别可变剪接时间的误差，比如：A3SS -> 剪接位点间至少要相差error个bp以上
     * @param silent 是否输出详细信息
     * @return list of SpliceEvent
     */
    fun identifyAS(error: Int, silent: Boolean = false): List<SpliceEvent> {

        if (!silent) {
            this.logger.info("Finding Alternative splicing events of ${this.chromosome}${this.strand}")
        }

        val res = mutableListOf<SpliceEvent>()

        val loggedEnds = HashSet<Sites>()

        // 识别可变剪接
        for ( i in this.starts.values.sorted() ) {
            if (i.size <= 1) {
                continue
            }

            for ( j in i.getSites() ) {
                val tmpEnds = this.ends[j.site]

                identifyBetweenSameStartEnd(
                        starts = i,
                        ends = tmpEnds,
                        res = res,
                        error = error
                )

                tmpEnds?.let {
                    loggedEnds.add(tmpEnds)
                }


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
                        res = res,
                        error = error
                )
            }
        }

        // finding MXE
        var sites = this.starts.keys.asSequence().sorted().toMutableList()
        sites.addAll(this.ends.keys)
        sites = sites.asSequence().distinct().sorted().toMutableList()

        val as35Starts = this.as35Sites["start"]!!.asSequence().sorted().distinct().toList()
        val as35Ends = this.as35Sites["end"]!!.asSequence().distinct().sorted().toList()

        var firstEnds = 0
        for ( i in 0..(as35Starts.size - 1)) {
            var j = firstEnds

            while ( j < as35Ends.size - 1 ) {

                // 如果j太小，就跳过
                if ( as35Ends[j] <= this.starts[as35Starts[i]]!!.getExtremeSite() ) {
                    if ( firstEnds != j ) {

                    }
                    j ++
                    firstEnds = j
                    continue
                }

                val currentStart = this.starts[as35Starts[i]]!!.getSites()
                val currentEnd = this.ends[as35Ends[j]]!!.getSites()

                for ( e1 in 0..(currentStart.size - 1) ) {
                    for ( e2 in (e1 + 1)..(currentStart.size - 1) ) {
                        for ( s1 in 0..(currentEnd.size - 1) ) {
                            for ( s2 in (s1 + 1)..(currentEnd.size - 1) ) {
                                val mxeSites = listOf(
                                        as35Starts[i],
                                        currentStart[e1].site,
                                        currentEnd[s1].site,
                                        currentStart[e2].site,
                                        currentEnd[s2].site,
                                        as35Ends[j]
                                )

                                if (
                                        this.isOrdered(mxeSites) &&
                                        sites.indexOf(currentStart[e2].site) - sites.indexOf(currentEnd[s1].site) == 1
                                ) {

                                    val tmpEvent = SpliceEvent(
                                            event = "MXE",
                                            chromosome = this.chromosome,
                                            start = as35Starts[i],
                                            end = as35Ends[j],
                                            strand = this.strand,
                                            sliceSites = mxeSites.toMutableList()
                                    )

                                    val values = when ( this.strand ) {
                                        '-' -> listOf(
                                                this.ends[as35Ends[j]]!!.getPsi( target = currentEnd[s1].site ),
                                                this.starts[as35Starts[i]]!!.getPsi( target = currentStart[e1].site )
                                        )
                                        else -> listOf(
                                                this.starts[as35Starts[i]]!!.getPsi( target = currentStart[e2].site ),
                                                this.ends[as35Ends[j]]!!.getPsi( target = currentEnd[s2].site )
                                        )
                                    }

                                    tmpEvent.psi = values.sum() / values.size
                                    res.add(tmpEvent)
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


    /**
     * 获取范围内所有符合要求的junctions位点，以start位点为标准
     * 获取范围内start对应的所有点的
     */
    fun getSites( range: Pair<Int, Int> ): List<Site> {
        val res = mutableListOf<Site>()

        for ( i in range.first..range.second ) {
            this.starts[i]?.let{
                res.addAll( it.getSites() )
            }
        }

        return res
    }


    override fun toString(): String {
        var output = ""
        for ( (start, value) in this.starts ) {
            for ( end in value.getSites()) {
                var tmp = "${this.chromosome}:$start-${end.site}${this.strand}\t${end.count}"

                if ( !end.source["gene"]!!.isEmpty() ) {
                    tmp = "$tmp\t${end.source["gene"]!!.joinToString(prefix = "", postfix = "", separator = ",")}"
                }

                if ( !end.source["transcript"]!!.isEmpty() ) {
                    tmp = "$tmp\t${end.source["transcript"]!!.joinToString(prefix = "", postfix = "", separator = ",")}"
                }

                output += "$tmp\n"
            }
        }
        return output
    }
}