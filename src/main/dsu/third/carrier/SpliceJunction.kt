package dsu.third.carrier

import dsu.carrier.Genes
import java.util.Objects

/**
 * @since 2018.06.21
 * @version 0.1
 * @author Zhang Yiming
 *
 * 可变剪接信息的载体
 */

class SpliceJunction(val gene: Genes) {

    data class Event(
            val name: String,
            val chromosome: String,
            val sites: List<Int>
    ) {
        /**
         * 重载toString
         * @return 名称，start和end以tab分割
         */
        override fun toString(): String {
            return "${this.name}\t${this.chromosome}:${this.sites.sorted().joinToString(separator = "-", prefix = "", postfix = "")}"
        }


        /**
         * 重载，识别两个事件是否相同
         */
        override fun equals(other: Any?): Boolean {
            try {
                return this.hashCode() == other!!.hashCode()
            } catch (err: NullPointerException) {
                return false
            }

        }


        /**
         * 重载，根据事件的信息，生成hashCode
         * @return hashCode
         */
        override fun hashCode(): Int {
            return Objects.hash(this.name, this.chromosome, this.sites.sorted())
        }

    }

    val events = mutableListOf<Event>()

    /**
     * 重载
     * @return 所有事件构成的字符串，可以直接写到文件，每行一个事件
     */
    override fun toString(): String {
        var output = ""

        for ( i in this.events.distinct() ) {
            output += "${this.gene.chromosome}:${this.gene.start}-${this.gene.end}\t${this.gene.geneId}\t${this.gene.transcriptId}\t$i\n"
        }

        return output
    }

    fun addEvent(name: String, chromosome: String, sites: List<Int>) {

        this.events.add(Event(
            name = name,
            chromosome = chromosome,
            sites = sites.sorted()
        ))
    }

    override fun hashCode(): Int {
        return Objects.hash(this.toString())
    }

    override fun equals(other: Any?): Boolean {
        if ( other == null ) {
            return false
        }
        return this.toString() == other.toString()
    }

}