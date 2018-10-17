package com.splice4k.base

import java.util.*


/**
 * @author Zhangyiming
 * @since 2018.09.05
 * @version 20180929
 * 二代算法记录相应的可变剪接的事件信息的类
 */


/**
 * @param event 事件类型
 * @param chromosome 染色体
 * @param start 起始位点
 * @param end 终止位点
 * @param strand 链
 * @param sliceSites 剪接位点
 */
class SpliceEvent(
        val event: String,
        chromosome: String,
        start: Int,
        end: Int,
        val strand: Char,
        sliceSites: MutableList<Int>
): GenomicLoci(chromosome, start, end) {

    val sliceSites = sliceSites
    get() {
        field.sort()
        return field
    }

    var psi: Double? = null
    val otherPSi = mutableListOf<Double>()
    var isNovel = true


    /**
     * 获取其他位点的PSI的值
     * @return string; None -> NA
     */
    fun getOtherPsi(): String {
        return when( this.otherPSi.isEmpty() ) {
            true -> "NA"
            false -> this.otherPSi.joinToString(separator = ",")
        }
    }


    override fun hashCode(): Int {
        return Objects.hash(this.event, this.chromosome, this.start, this.end, this.sliceSites.sorted())
    }

    override fun equals(other: Any?): Boolean {
        return this.hashCode() == other?.hashCode()
    }

    override fun toString(): String {
        return "${this.chromosome}:${this.sliceSites.first()}-${this.sliceSites.last()}${this.strand}\t${this.event}\t${this.chromosome}:${this.sliceSites.joinToString(prefix = "", postfix = "", separator = "-")}"
    }

}
