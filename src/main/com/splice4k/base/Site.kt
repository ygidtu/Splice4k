package com.splice4k.base


/**
 * @since 2018.09.21
 * @version 20180921
 * @author Zhang yiming
 *
 * 记录单个位点的位置，来源和数量信息
 */


/**
 * @param site 位点
 */
class Site(val site: Int): Comparable<Site> {
    var source = hashMapOf(
            "transcript" to mutableListOf<String>(),
            "gene" to mutableListOf()
    )
    var count = 0

    /**
     * 添加位点来源
     * @param gene 基因id
     * @param transcript 转录本id
     */
    fun addSource( gene: String?=null, transcript: String? = null ) {
        if ( gene != null ) {
            this.source["gene"]!!.add(gene)
        }

        if ( transcript != null ) {
            this.source["transcript"]!!.add(transcript)
        }
    }

    /**
     * 添加位点数量
     * @param count 频数，默认1
     */
    fun addCount(count: Int? = 1) {

        this.count += when(count) {
            null -> 1
            else -> count
        }
    }


    override fun hashCode(): Int {
        return this.site
    }


    override fun equals(other: Any?): Boolean {
        return this.hashCode() == other?.hashCode()
    }


    override fun compareTo(other: Site): Int {
        return this.site - other.site
    }

    override fun toString(): String {
        return "${this.site}[${this.count}]"
    }
}

