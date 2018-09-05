package dsu.carrier

import java.util.Objects

/**
 * @author zhangyiming
 * @since 2018.06.20
 * @version 0.1
 * 基因与Reads的配对
 */

class GeneRead(val gene: Genes, val reads: Genes): Comparable<GeneRead> {

    var overlap: Int = this.gene.getOverlap(this.reads) ?: 0

    // 比例改为相对于较短的那一个来算
    var overlapPercent: Double = this.gene.overlapPercent(this.reads)


    /**
     * 判断基因和reads的exon的重合程度是否符合90%的要求
     * 任意一对符合要求即可
     * @return Boolean  true符合；false不符合
     */
    fun isGeneReadsExonsOverlapQualified(overlapStandard: Double = 90.0): Boolean {
        var i = 0; var j = 0
        val geneExons = this.gene.exons.sorted()
        val readsExons = this.reads.exons.sorted()

        while (i < geneExons.size && j < readsExons.size ) {
            val tmpGene = geneExons[i]
            val tmpRead = readsExons[j]

            when {
                tmpGene.isUpStream(tmpRead) -> i++
                tmpGene.isDownStream(tmpRead) -> j++
                else -> {
                    if (tmpGene.overlapPercent(tmpRead, true) >= overlapStandard) {
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
        return Objects.hash(this.gene.transcriptId, this.reads.transcriptId)
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


    /**
     * toString是按照gene|reads的顺序构建
     * 这个function反过来
     * @return string
     */
    fun toStringReverse() : String {
        if (this.gene.transcriptId == "." && this.gene.geneName == "." && this.gene.length <= 0 ) {
            return "${this.reads}|None"
        }
        return "${this.reads}|${this.gene}"
    }
}