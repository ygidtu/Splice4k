package dsu.carrier

import java.util.Objects


/**
 * @author Zhangyiming
 * @since 2018.09.05
 * @version 20180905
 * 二代算法记录相应的可变剪接的事件信息的类
 */

class SpliceEvent(
        val event: String,
        chromosome: String,
        start: Int,
        end: Int,
        val strand: Char,
        sliceSites: MutableList<Int>
): GenomicLoci(chromosome, start, end) {

    val sliceSites = sliceSites
    get() {
        field.sort()
        return field
    }


    override fun hashCode(): Int {
        return Objects.hash(this.event, this.chromosome, this.start, this.end, this.sliceSites.sorted())
    }

    override fun equals(other: Any?): Boolean {
        return this.hashCode() == other?.hashCode()
    }

    override fun toString(): String {
        return "${this.chromosome}:${this.sliceSites.first()}-${this.sliceSites.last()}${this.strand}\t${this.event}\t${this.chromosome}:${this.sliceSites.joinToString(prefix = "", postfix = "", separator = "-")}"
    }

}
