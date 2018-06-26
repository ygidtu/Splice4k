package main.carrier

import java.util.Objects

/**
 * @author zhangyiming
 * @since 2018.06.20
 * @version 0.1
 * 基因与Reads的配对
 */

class Template(
        val gene: Genes,
        val template: Genes
): Comparable<Template> {
    var geneExons: MutableList<Array<Int>>

    init {
        this.geneExons = this.filterExons()
    }


    /**
     * 检查两个外显子是否有重合位点
     * @param first 外显子位点 Array(start, end)
     * @param second 外显子位点 Array(start, end)
     * @return true 有。false 没有
     */
    private fun isExonOverlap(first: Array<Int>, second: Array<Int>): Boolean {
        return (first[0] <= second[1] && first[1] >= second[0])
    }


    /**
     * 检查第一个外显子是否在第二个的上游
     * @param first 外显子位点 Array(start, end)
     * @param second 外显子位点 Array(start, end)
     * @return true 有。false 没有
     */
    private fun isExonUpStream(first: Array<Int>, second: Array<Int>): Boolean {
        return (first[1] < second[0])
    }


    /**
     * 检查第一个外显子是否在第二个的下游
     * @param first 外显子位点 Array(start, end)
     * @param second 外显子位点 Array(start, end)
     * @return true 有。false 没有
     */
    private fun isExonDownStream(first: Array<Int>, second: Array<Int>): Boolean {
        return (first[0] > second[1])
    }


    /**
     * 过滤一遍基因的外显子
     * 仅保留基因中与template中有重合的外显子
     * @return 一列外显子的坐标
     */
    private fun filterExons(): MutableList<Array<Int>> {
        val geneExons = mutableListOf<Array<Int>>()

        var i = 0
        var j = 0

        while (i < this.gene.exons.size && j < this.template.exons.size) {
            val tmpGene = this.gene.exons[i]
            val tmpRead = this.template.exons[j]


            if (this.isExonUpStream(tmpGene, tmpRead)) {
                j++
            } else if (this.isExonDownStream(tmpGene, tmpRead)) {
                i++
            } else {
                geneExons.add(tmpGene)
                j++
            }
        }

        return geneExons
    }

    /**
     * toString 重载，将对应的基因和template配对成一行文本输出
     * @return gene|read
     */
    override fun toString(): String {
        var exonString = ""
        for (i in this.template.exons) {
            exonString += i.joinToString(prefix = "", postfix = ",", separator = ",")
        }

        return "${this.template.chrom}\t${this.template.start}\t${this.template.end}\t${this.gene.transcriptId}\t${this.gene.geneName}\t$exonString"
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