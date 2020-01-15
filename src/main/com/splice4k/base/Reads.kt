package com.splice4k.base

import com.splice4k.base.Exons
import kotlin.math.max
import kotlin.math.min
import java.util.Objects

// 数据类，用以提取BAM，gmap align文件中的信息，构建AS所需图
class Reads(
        chromosome: String,
        start: Int,
        end: Int,
        val strand: Char,
        val name: String,
        val introns: List<Int>
): GenomicLoci(chromosome, start, end) {


    override fun toString(): String {
        return "${this.chromosome}:${this.start}-${this.end}:${this.strand}\t${this.name}"
    }

    fun getExon(): List<Exons> {
        val res = mutableListOf<Exons>()

        for (i in 0.until(this.introns.size)) {

            var start = -1
            var end = -1

            if (i == 0) {
                start = this.start
                end = this.introns[0] - 1
            } else if (i == this.introns.size - 1) {
                start = this.introns.last() + 1
                end = this.end
            } else if (i % 2 == 1) {
                start = this.introns[i]
                end = this.introns[i + 1]
            }

            if (start > 0 && end > 0) {
                res.add(Exons(
                        chromosome=this.chromosome,
                        start = start,
                        end = end,
                        strand = this.strand,
                        exonId = "${this.name}.${i / 2 + 1}"
                ))
            }
        }

        return res
    }

    /**
     * 重载hashCode
     */
    override fun hashCode(): Int {
        return Objects.hash(this.chromosome, this.start, this.end, this.strand, this.introns)
    }

    /**
     * 重载equals
     */
    override fun equals(other: Any?): Boolean {
        return other is Reads && this.hashCode() == other.hashCode()
    }


    /**
     * 重载的compareTo，用来比较位点的上下游
     */
    override fun compareTo(other: GenomicLoci): Int {

        return if (other !is Reads) {
            compareValuesBy(this, other, {it.chromosome}, {it.start}, {it.end})
        } else if (compareValuesBy(this, other, {it.chromosome}, {it.start}, {it.end}, {it.getExon().size}) != 0) {
            compareValuesBy(this, other, {it.chromosome}, {it.start}, {it.end}, {it.getExon().size})
        } else {
            var res = 0
            val thisExons = this.getExon()
            val otherExons = other.getExon()

            for (i in 0.until(thisExons.size - 1)) {
                if (thisExons[i] != otherExons[i]) {
                    res = thisExons[i].compareTo(otherExons[i])
                    break
                }
            }
            res
        }
    }

    fun addReads(other: Reads) {
        this.start = min(this.start, other.start)
        this.end = max(this.end, other.end)
    }
}


