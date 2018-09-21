package dsu.carrier

import java.util.*

/**
 * @author zhangyiming
 * @since 2018.07.04
 * @version 20180903
 * 外显子常用的function 三连
 */



class Exons( start: Int, end: Int): GenomicLoci(start, end) {
    var source: String = ""
    var strand: Char = '.'


    /**
     * 构造器
     * @param chromosome 染色体
     * @param start 起始位点
     * @param end 终止位点
     * @param source 来源，比如transcript Id
     */
    constructor(
            chromosome: String,
            start: Int,
            end: Int,
            source: String
    ): this(start, end) {
        this.chromosome = chromosome
        this.source = source
    }


    /**
     * 构造器
     * @param chromosome 染色体
     * @param start 起始位点
     * @param end 终止位点
     * @param strand 正负链
     */
    constructor(
            chromosome: String,
            start: Int,
            end: Int,
            strand: Char
    ): this(start, end) {
        this.chromosome = chromosome
        this.strand = strand
    }


    /**
     * 构造器
     * @param chromosome 染色体
     * @param start 起始位点
     * @param end 终止位点
     * @param source 来源，比如transcript Id
     * @param strand 正负链
     */
    constructor(
            chromosome: String,
            start: Int,
            end: Int,
            source: String,
            strand: Char
    ): this(chromosome, start, end, source) {
        this.strand = strand
    }


    /**
     * 获取某个exon 起始位点或者终止位点的hash，便于构建same-start和same-end
     * @param start 是否返回start的hash：true -> start的；false -> end
     * @return hashCode of chromosome, one of start or end, strand
     */
    fun getSiteHash(start: Boolean = true): Int {
        return when(start) {
            true -> Objects.hash(this.chromosome, this.start, this.strand)
            else -> Objects.hash(this.chromosome, this.end, this.strand)
        }
    }
}




