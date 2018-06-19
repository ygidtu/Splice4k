package informationCarrier

import jdk.nashorn.internal.runtime.regexp.joni.exception.ValueException

/**
 * @since 2018.06.14
 * @author zhangyiming
 * @version 0.1
 * 基因信息载体
*/

class Genes: Comparable<Genes> {

    var fileType : String

    var chrom: String = "."

    var start: Int = -1
    set(value) {
        if (value > 0) field = value
        else throw ValueException("${value} < 0")
    }

    var end: Int = -1
    set(value) {
        if (value >= this.start && value > 0) field = value
        else throw ValueException("End (${value}) < start ${this.start} or ${value} < 0")
    }

    var strand: Char = '.'
    set(value) {
        if (value in arrayOf('+', '-', '.')) field = value
    }

    var geneId: String = "."
    var geneName: String = "."

    var transcriptId: String? = null

    var exons: MutableList<Array<Int>> = mutableListOf()

    /**
     * constructors
     * @chrom chromosome
     * @start start site
     * @end end site
     * @geneName gene name
     * @information gtf和gff文件中额外的信息封装成的map
     * @strand + -
     */
    constructor(chrom: String, start: Int, end: Int) {
        this.chrom = chrom
        this.start = start
        this.end = end
    }

    constructor(chrom: String, start: Int, end: Int, geneName: String) : this(chrom, start, end) {
        this.geneName = geneName
    }

    constructor(chrom: String, start: Int, end: Int, information: Map<String, String>) : this(chrom, start, end) {
        for ((k, v) in information) {
            when (k) {
                "gene_id" -> this.geneId = v
                "ID" -> this.geneId = v
                "gene_name" -> this.geneName = v
                "Name" -> this.geneName = v
                "transcript_id" -> this.transcriptId = v
            }
        }
    }

    constructor(chrom: String, start: Int, end: Int, strand: Char, information: Map<String, String>) :
            this(chrom, start, end, information) {
        this.strand = strand
    }

    override fun toString(): String {
        return "${this.chrom}\t${this.start}\t${this.end}\t${this.geneId}\t${this.geneName}\t${this.strand}"
    }

    fun isGtf(): Boolean {
        if (this.fileType == "gtf") {
            return true
        }
        return false
    }

    override fun equals(other: Any?): Boolean {
        return this.transcriptId == other.toString()
    }

    /**
     * 首选transcript id作为hashCode的来源，其次选择gene id
     */
    override fun hashCode(): Int {
        return this.transcriptId?.hashCode() ?: this.geneId.hashCode()
    }

    /**
     * 重载compareTo函数，用来比较两个位点之间是否存在对应的重合
     * 也便于进行排序
     * @other 另外一个Genes类
     */
    override fun compareTo(other: Genes): Int {
        if (this.chrom > other.chrom) {
            return 1
        } else if (this.chrom < other.chrom) {
            return -1
        } else {
            if (this.start > other.start) {
                return 1
            } else if (this.start < other.start) {
                return -1
            } else {
                if (this.end > other.end) {
                    return 1
                } else if (this.end < other.end) {
                    return -1
                } else {
                    return 0
                }
            }
        }


    }

    /**
     * 添加exon
     */
    fun addExons(exon: Array<Int>) {
        this.exons.add(exon.subtract(2))
    }
}