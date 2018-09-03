package dsu.carrier


import java.util.*



/**
 * @since 2018.06.14
 * @author zhangyiming
 * @version 20180903
 * 基因信息载体
*/

class Genes(chromosome: String, start: Int, end: Int): GenomicLoci(chromosome, start, end) {

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
                '.' -> "${this.chromosome}:${this.start}-${this.end}"
                else -> "${this.chromosome}:${this.start}-${this.end}${this.strand}"
            }
            else -> {
                return field
            }
        }
    }

    var exons: MutableList<Exons> = mutableListOf()
    get() {
        return field.sorted().toMutableList()
    }

    var parent: String = ""
    get() {
        return when (field) {
            "" -> this.geneId
            else -> field
        }
    }

    val length: Int
    get() {
        return this.end - this.start
    }


    /**
     * constructors
     * @chromosome 染色体名称
     * @start 起始位点
     * @end 终止位点
     * @geneName gene name
     * @information gtf和gff文件中额外的信息封装成的map
     * @strand + -
     */
    constructor(
            chromosome: String,
            start: Int,
            end: Int,
            geneName: String
    ): this(chromosome, start, end) {
        this.geneName = geneName
    }

    /**
     * constructors
     * @chromosome 染色体名称
     * @start 起始位点
     * @end 终止位点
     * @geneName gene name
     * @strand + -
     */
    constructor(
            chromosome: String,
            start: Int,
            end: Int,
            geneName: String,
            strand: Char
    ) :
            this(chromosome, start, end, geneName) {
        this.strand = strand
    }

    /**
     * constructors
     * @chromosome 染色体名称
     * @start 起始位点
     * @end 终止位点
     * @information gtf和gff文件中额外的信息封装成的map
     */
    constructor(
            chromosome: String,
            start: Int,
            end: Int,
            information: Map<String, String>
    ) : this(chromosome, start, end) {
        for ((k, v) in information) {
            when (k) {
                "gene_id" -> this.geneId = v
                "Parent" -> this.geneId = v
                "gene_name" -> this.geneName = v
                "Name" -> this.geneName = v
                "transcript_id" -> this.transcriptId = v
                "ID" -> this.transcriptId = v
            }
        }
    }

    /**
     * constructors
     * @chromosome 染色体名称
     * @start 起始位点
     * @end 终止位点
     * @information gtf和gff文件中额外的信息封装成的map
     * @strand + -
     */
    constructor(chromosome: String, start: Int, end: Int, strand: Char, information: Map<String, String>) :
            this(chromosome, start, end, information) {
        this.strand = strand
    }

    override fun equals(other: Any?): Boolean {
        return this.transcriptId == other.toString()
    }

    /**
     * 首选transcript id作为hashCode的来源，其次选择gene id
     */
    override fun hashCode(): Int {
        return Objects.hash(this.transcriptId)
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
    override fun compareTo(other: GenomicLoci): Int {
        return when (other) {
            is Genes ->  compareValuesBy( this, other, {it.chromosome}, {it.start}, {it.end}, {it.strand})
            else -> compareValuesBy(this, other, {it.chromosome}, {it.start}, {it.end})
        }
    }

    /**
     * 返回目前已有的数据
     * @return list of message
     */
    fun get() : List<String> {
        var exonString = ""
        for (i in this.exons) {
            exonString += "$i,"
        }

        return listOf(
                this.chromosome,
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

/*
fun main(args: Array<String>) {
    val test = Genes(
            "chr1", 1, 100,
            "test"
    )

    test.addExons(arrayOf(10, 20, 30))

    println(test)
}
*/