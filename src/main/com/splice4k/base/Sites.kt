package com.splice4k.base


/**
 * @since 2018.09.21
 * @version 20180921
 * @author Zhang yiming
 * 内部类，记录所有位点信息
 */

/**
 * @param node same start或者same end的核心（就是same的那个start或者end）
 */
class Sites( val node: Int ): Comparable<Sites> {
    private val pos = mutableMapOf<Int, Site>()

    var total: Double = 0.toDouble()

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

        this.total += freq ?: 1

        this.pos[site] = tmp
    }

    /**
     * 按照index获取Sites中的位点
     * @param target 获取某个index位置上的Site，输入的index
     * @return Site？
     */
    fun getSiteByIndex(target: Int): Site {
        return this.pos[this.pos.keys.sorted()[target]]!!
    }

    /**
     * 获取sites中记录的所有位点
     * @return 获取这个节点上对应的所有其他点
     */
    fun getSites(): List<Site> {
        return this.pos.values.toMutableList().sorted()
    }


    /**
     * 获取某个特定位点的PSI值
     * @param target 所需PSI值的位点
     * @return Double? null -> means site doesn't exists; double -> PSI value
     */
    fun getPsi(target: Int): Double {
        this.pos[target]!!.let {
            return it.count / this.total
        }
    }


    /**
     * retrive count value of specific site
     * @param target the site
     * @return Int -> count of that site
     */
    fun getCount(target: Int): Int {
        return this.pos[target]!!.count
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

}