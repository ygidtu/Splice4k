package com.splice4k.base


/**
 * @author zhangyiming
 * @since 2018.07.04
 * @version 20180903
 * 外显子常用的function 三连
 */


/**
 * @param chromosome 染色体
 * @param start 起始位点
 * @param end 终止位点
 * @param exonId 外显子id
 */
class Exons(
        chromosome: String,
        start: Int,
        end: Int,
        val exonId: String

): GenomicLoci(chromosome, start, end) {

    var source = mutableMapOf(
            "gene" to mutableSetOf<String>(),
            "transcript" to mutableSetOf()
    )

    // 专门用于isoform组装，表明来源所用
    var annotation = ""

    override fun toString(): String {
        return this.source["gene"]!!.asSequence().distinct().joinToString(prefix = "", postfix = "\t", separator = ",") +
                this.source["transcript"]!!.asSequence().distinct().joinToString(prefix = "", postfix = "\t", separator = ",") +
                this.exonId
    }
}




