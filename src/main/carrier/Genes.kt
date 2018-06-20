package main.carrier

import jdk.nashorn.internal.runtime.regexp.joni.exception.ValueException

/**
 * @since 2018.06.14
 * @author zhangyiming
 * @version 0.1
 * 基因信息载体
*/

class Genes: Comparable<Genes> {

    var fileType : String = "gtf"

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
    set(value) {
        field = when(":" in value) {
            true -> value.split(":")[1]
            false -> value
        }
    }

    var geneName: String = "."

    var transcriptId: String = "."
    get() {
        return when (field) {
            "." -> when(this.strand) {
                '.' -> "${this.chrom}:${this.start}-${this.end}"
                else -> "${this.chrom}:${this.start}-${this.end}${this.strand}"
            }
            else -> this.geneId
        }
    }

    var exons: MutableList<Array<Int>> = mutableListOf()
    get() {
        return field.sortedWith(compareBy({it[0]})).toMutableList()
    }

    var parent: String = ""


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

    constructor(chrom: String, start: Int, end: Int, geneName: String, strand: Char) :
            this(chrom, start, end, geneName) {
        this.strand = strand
    }

    constructor(chrom: String, start: Int, end: Int, information: Map<String, String>) : this(chrom, start, end) {
        for ((k, v) in information) {
            when (k) {
                "gene_id" -> this.geneId = v
                "ID" -> this.geneId = v
                "gene_name" -> this.geneName = v
                "Name" -> this.geneName = v
                "transcript_id" -> this.transcriptId = v
                "Parent" -> this.parent = v
            }
        }
    }

    constructor(chrom: String, start: Int, end: Int, strand: Char, information: Map<String, String>) :
            this(chrom, start, end, information) {
        this.strand = strand
    }

    constructor(chrom: String, start: Int, end: Int, geneName: String, information: Map<String, String>) :
            this(chrom, start, end, information) {
        this.geneName = geneName
    }

    constructor(chrom: String, start: Int, end: Int, geneName: String, strand: Char, information: Map<String, String>) :
            this(chrom, start, end, geneName, information) {
        this.strand = strand
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
        return this.transcriptId.hashCode()
    }

    /**
     * toString
     */
    override fun toString(): String {
        return this.get().joinToString(separator = "\t", prefix = "", postfix = "")
    }

    /**
     * 重载compareTo函数，用来比较两个位点之间是否存在对应的重合
     * 也便于进行排序
     * @other 另外一个Genes类
     */
    override fun compareTo(other: Genes) = compareValuesBy(this, other, {it.chrom}, {it.start}, {it.end})

    /**
     * 添加exon
     * @param exon 两个元素的数据，int，分别为外显子的起始和终止位点
     */
    fun addExons(exon: Array<Int>) {
        this.exons.add(exon.sliceArray(0..1))
    }

    /**
     * 判断这个位点是否覆盖另外一个位点
     * @other Genes
     */
    fun isCover(other: Genes): Boolean {
        if (
                this.chrom == other.chrom &&
                this.start <= other.start &&
                this.end >= other.end
        ) {
            return true
        }
        return false
    }

    /**
     * 判断该位点是否为另一个位点的上游
     * @other Genes
     */
    fun isUpStream(other: Genes) : Boolean {
        if (this.chrom != other.chrom) {
            return this.chrom < other.chrom
        }

        return this.end < other.start
    }


    /**
     * 判断该位点是否在另一个位点的下游
     * @other Genes
     */
    fun isDownStream(other: Genes) : Boolean {
        if (this.chrom != other.chrom) {
            return this.chrom > other.chrom
        }

       return this.start > other.end
    }

    /**
     * 返回目前已有的数据
     * @return list of message
     */
    fun get() : List<String> {
        var exonString = ""
        for (i in this.exons) {
            exonString += i.joinToString(prefix = "", postfix = "\t", separator = "\t")
        }

        return listOf(
                this.chrom,
                this.start.toString(),
                this.end.toString(),
                this.strand.toString(),
                this.geneId,
                this.geneName,
                this.transcriptId,
                exonString
        )
    }
}


fun main(args: Array<String>) {
    val test = Genes(
            "chr1", 1, 100,
            "test"
    )

    test.addExons(arrayOf(10, 20, 30))

    println(test)
}