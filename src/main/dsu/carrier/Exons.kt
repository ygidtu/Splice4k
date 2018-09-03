package dsu.carrier

import kotlin.math.*
/**
 * @author zhangyiming
 * @since 2018.07.04
 * @version 20180903
 * 外显子常用的function 三连
 */



class Exons( start: Int, end: Int): GenomicLoci(start, end) {
    var source: String = ""
    var strand: Char = '.'

    constructor(
        chromosome: String,
        start: Int,
        end: Int,
        source: String
    ): this(start, end) {
        this.chromosome = chromosome
        this.source = source
    }

    constructor(
        chromosome: String,
        start: Int,
        end: Int,
        source: String,
        strand: Char
    ): this(chromosome, start, end, source) {
        this.strand = strand
    }
}




