package splice4k.smrt.base

import splice4k.base.Genes
import java.util.*

/**
 * @author Zhang yiming
 * @since 2018.06.20
 * @version 20180926
 * 基因与Reads的配对
 */

class GeneRead(val gene: Genes, val reads: Genes): Comparable<GeneRead> {

    // 比例改为相对于较短的那一个来算
    var overlapPercent: Double = this.gene.overlapPercent(this.reads)


    /**
     * 判断基因和reads的exon的重合程度是否符合90%的要求
     * 任意一对符合要求即可
     * @return Boolean  true符合；false不符合
     */
    fun isGeneReadsExonsOverlapQualified(error: Int = 3): Boolean {
        val geneExons = this.gene.exons
        val readsExons = this.reads.exons

        var i = 0; var j = 0; var match = 0
        while (i < geneExons.size && j < readsExons.size) {
            when  {
                geneExons[i] < readsExons[j] - error -> i ++
                geneExons[i] > readsExons[j] + error -> j ++
                else -> {
                    match ++

                    if (match >= 2) {
                        return true
                    }
                    j++
                }
            }
        }
        return false
    }

    /**
     * toString 重载，将对应的基因和reads配对成一行文本输出
     * @return gene|read
     */
    override fun toString(): String {
        if (this.gene.transcriptId == "." && this.gene.geneName == "." && this.gene.length <= 0 ) {
            return "None|${this.reads}"
        }
        return "${this.gene}|${this.reads}"
    }

    /**
     * hashCode重载
     */
    override fun hashCode(): Int {
        return Objects.hash(this.gene.toString(), this.reads.toString())
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
     * 先根据基因的id号进行比较，再根据reads的位点进行排序。
     * 保证同一个基因的所有配对结果放在一起，然后进行比较的时候又省掉排序的过程
     * @param other 另外一个GeneRead类
     * @return -1 小于, 0等于, 1大于
     */
    override fun compareTo(other: GeneRead): Int {
        return when {
            this.gene > other.gene -> 1
            this.gene < other.gene -> -1
            else -> {
                when {
                    this.reads > other.reads -> 1
                    this.reads < other.reads -> -1
                    else -> 0
                }
            }
        }
    }
}