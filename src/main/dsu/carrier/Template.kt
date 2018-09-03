package dsu.carrier


import java.util.Objects

/**
 * @author zhangyiming
 * @since 2018.06.20
 * @version 20180903
 * 基因与Reads的配对
 */

class Template(
        val gene: Genes,
        val template: Genes
): Comparable<Template> {
    var geneExons = this.gene.exons.sorted().toMutableList()
    get() {
        return field.sorted().toMutableList()
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

        return "${this.template.chromosome}\t${this.template.start}\t${this.template.end}\t${this.gene.transcriptId}\t${this.gene.geneName}\t$exonString"
    }

    /**
     * hashCode重载
     */
    override fun hashCode(): Int {
        return Objects.hash(this.gene.transcriptId, this.template.transcriptId)
    }

    /**
     * equals重载
     */
    override fun equals(other: Any?): Boolean {
        if (other != null) {
            return this.hashCode() == other.hashCode()
        }
        return false
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
            this.gene > other.gene -> 1
            this.gene < other.gene -> -1
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