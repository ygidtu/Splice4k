package main.Carrier

import java.util.Objects

/**
 * @author zhangyiming
 * @since 2018.06.20
 * @version 0.1
 * 基因与Reads的配对
 */

class GeneRead(val gene: Genes, val reads: Genes) {

    /**
     * toString 重载，将对应的基因和reads配对成一行文本输出
     * @return gene|read
     */
    override fun toString(): String {
        return "${this.gene}|${this.reads}"
    }

    /**
     *
     */
    override fun hashCode(): Int {
        return Objects.hash(this.gene.transcriptId, this.reads.transcriptId)
    }
}