package dsu.carrier


/**
 * @author zhangyiming
 * @since 2018.07.04
 * @version 20180903
 * 外显子常用的function 三连
 */



class Exons(

        chromosome: String,
        start: Int,
        end: Int,
        val strand: Char,
        val exonId: String

): GenomicLoci(chromosome, start, end) {

    var source = mutableMapOf(
            "gene" to "",
            "transcript" to ""
    )
}




