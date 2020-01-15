package com.splice4k.tools

import com.splice4k.base.Reads
import com.splice4k.base.GenomicLoci
import kotlin.math.abs

class RelativePos {



    fun isUpStream(first: GenomicLoci, second: GenomicLoci, error:Int = 3): Boolean {
        if (first.chromosome != second.chromosome) {
            return first.chromosome < second.chromosome
        }

        if (abs(first.start - second.start) >= error) {
            return first.start < second.start
        }

        if (abs(first.end - second.end) >= error) {
            return first.end < second.end
        }

        return false
    }


    fun isDownStream(first: GenomicLoci, second: GenomicLoci, error:Int = 3): Boolean {
        if (first.chromosome != second.chromosome) {
            return first.chromosome > second.chromosome
        }

        if (abs(first.start - second.start) >= error) {
            return first.start > second.start
        }

        if (abs(first.end - second.end) >= error) {
            return first.end > second.end
        }

        return false
    }


    fun isSame(first: GenomicLoci, second: GenomicLoci, error:Int = 3): Boolean {
        return first.chromosome == second.chromosome &&
                abs(first.start - second.start) < error &&
                abs(first.end - second.end) < error
    }

    fun isSameReads(first: Reads, second: Reads, error: Int = 3): Boolean {
        if(first.chromosome == second.chromosome &&
                abs(first.start - second.start) < error &&
                abs(first.end - second.end) < error &&
                first.strand == second.strand &&
                first.getExon().size == second.getExon().size) {

            val firstExons = first.getExon()
            val secondExons = second.getExon()

            for(i in 0.until(firstExons.size - 1)) {
                if (isSame(firstExons[i], secondExons[i], error)) {
                    return true
                }
            }
        }

        return false
    }

}