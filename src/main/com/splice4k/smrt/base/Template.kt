package splice4k.smrt.base

import splice4k.base.Genes
import java.util.*


/**
 * @author zhang yiming
 * @since 2018.06.20
 * @version 20180926
 * 基因与Reads的配对
 */

class Template(
        val template: Genes,
        reads: MutableList<Genes>
): Comparable<Template> {
    val reads = reads
    get() {
        field.sort()
        return field
    }

    /**
     * 获取这个组装好的模板对应的所有reads的exons
     * @return sorted list of exons
     */
    fun getReadsExons(): List<List<Int>> {
        val res = mutableListOf<List<Int>>()

        for ( i in this.reads.distinct() ) {
            res.add(i.exons)
        }

        return res
    }


    /**
     * toString 重载，将对应的基因和template配对成一行文本输出
     * @return gene|read
     */
    override fun toString(): String {
        var exonString = ""
        for (i in this.template.exons) {
            exonString += "$i,"
        }

        return "${this.template.chromosome}\t${this.template.start}\t${this.template.end}\t${this.template.transcriptId}\t${this.template.geneName}\t$exonString"
    }

    /**
     * hashCode重载
     */
    override fun hashCode(): Int {
        return Objects.hash(this.template.transcriptId, this.template.transcriptId)
    }

    /**
     * equals重载
     */
    override fun equals(other: Any?): Boolean {
        return this.hashCode() == other?.hashCode()
    }

    /**
     * 重载，用于作比较
     * 先根据基因的id号进行比较，再根据template的位点进行排序。
     * 保证同一个基因的所有配对结果放在一起，然后进行比较的时候又省掉排序的过程
     * @param other 另外一个GeneRead类
     * @return -1 小于, 0等于, 1大于
     */
    override fun compareTo(other: Template): Int {
        return when {
            this.template > other.template -> 1
            this.template < other.template -> -1
            else -> {
                when {
                    this.template > other.template -> 1
                    this.template < other.template -> -1
                    else -> 0
                }
            }
        }
    }
}