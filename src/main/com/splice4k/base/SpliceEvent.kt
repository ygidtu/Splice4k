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
    var isNovel = true
    var subtypes = "NA"


    /**
     * get PSI value in string format
     * only keep 3 decimals
     * @return PSI value in String, null -> NA; else -> 3 decimals
     */
    fun getPsi(): String {
        return when (this.psi) {
            null -> "NA"
            else -> String.format("%.3f", this.psi)
        }
    }

    /**
     * calculate hash code only based on the event type, genomic loci and splice sites
     * @return hash code
     */
    override fun hashCode(): Int {
        return Objects.hash(this.event, this.chromosome, this.start, this.end, this.sliceSites.sorted())
    }


    /**
     * does two splice events equal?
     * @return bool
     */
    override fun equals(other: Any?): Boolean {
        return this.hashCode() == other?.hashCode()
    }


    /**
     * format the spliced region, event type, subtype and splice sites into String
     * @return String, eg: chr1:100-300+\tA3\tNA\tchr1:100-100-200-300
     */
    override fun toString(): String {
        return "${this.chromosome}:${this.sliceSites.first()}-${this.sliceSites.last()}${this.strand}" +
                "\t${this.event}" +
                "\t${this.subtypes}" +
                "\t${this.chromosome}:${this.sliceSites.joinToString(separator = "-")}"
    }


    /**
     * get formatted string, just add new subtypes
     */
    fun getString(subtypes: String): String {
        return "${this.chromosome}:${this.sliceSites.first()}-${this.sliceSites.last()}${this.strand}" +
                "\t${this.event}" +
                "\t$subtypes" +
                "\t${this.chromosome}:${this.sliceSites.joinToString(separator = "-")}"
    }

}
