package dsu.carrier

import java.lang.NullPointerException

/**
 * @since 2018.09.21
 * @version 20180921
 * @author Zhang yiming
 * 内部类，记录所有位点信息
 */
class Sites( val node: Int ): Comparable<Sites>, Iterator<Int> {
    val pos = mutableMapOf<Int, Site>()

    var current = 0
    var size = this.pos.size
        get() {
            return this.pos.size
        }

    /**
     * 添加同start或者同end的位点
     * @param site int 位点
     * @param freq int? 位点出现的次数
     * @param transcript 位点来源的转录本
     * @param gene 位点来源的基因
     */
    fun addSite( site: Int, freq: Int? = null, transcript: String? = null, gene: String? = null ) {
        val tmp = when (this.pos.containsKey(site)) {
            true -> this.pos[site]!!
            else -> Site(site)
        }
        tmp.addSource(transcript = transcript, gene = gene)
        tmp.addCount(freq)

        this.pos[site] = tmp
    }


    fun getSites(): List<Site> {
        return this.pos.values.toList()
    }


    /**
     * 获取离node最近的点
     */
    fun getClosestSite(): Int {
        val tmp = this.pos.values.sorted()

        return when(
            kotlin.math.abs(tmp.first().site - this.node) <
                    kotlin.math.abs(tmp.last().site - this.node)
            ) {
            true -> tmp.first().site
            else -> tmp.last().site
        }
    }


    /**
     * 获取最长的那个位点，即距离node最远的位点
     */
    fun getExtremeSite(): Int {
        val tmp = this.pos.values.sorted()

        return when(
            kotlin.math.abs(tmp.first().site - this.node) >
                    kotlin.math.abs(tmp.last().site - this.node)
            ) {
            true -> tmp.first().site
            else -> tmp.last().site
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
        return this.pos.values.sorted()[current - 1].site
    }

    override fun hasNext(): Boolean {
        return this.current < this.pos.size
    }

}