package com.splice4k.tools

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File

/**
 * 专门用于IR事件中PSI的计算
 * @author Zhang Yiming
 * @since 2018.09.29
 * @version 20180929
 */



class PsiOfIR() {
    /**
     * private function
     * 从单条bam中提取cigar中的N区域
     * @param record 单条SAM/BAM的信息
     * @return 列表，记录了所有的intron的边界信息
     */
    private fun extractSpliceFromCigar( record: SAMRecord ): List<Int> {
        var position = record.alignmentStart
        val results = mutableListOf( position )
        val tmp = mutableListOf<Char>()

        for (i in record.cigar.toString()) {
            if (i in '0'..'9') {  // 如果是数字，就加到list中
                tmp.add(i)
            } else {
                if (tmp.size == 0) {
                    continue
                }

                position += tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()


                if ( i in arrayOf('D', 'H', 'N') ) {
                    results.add(
                            position - tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                    )

                    results.add(position - 1)
                }
                tmp.clear()
            }
        }

        results.add( record.alignmentEnd )
        return results
    }


    /**
     * 计算IR的PSI值
     * 该值与其他的PSI值不同，该值需要的是
     * - 前后两个外显子范围20bp范围内，每个bp的coverage的平均值
     * - junctions区域平均的coverage
     * 然后junctions / ((exon1 + exons2) / 2)
     */
    fun getPsi(
            chromosome: String,
            regionStart: Int,
            regionEnd: Int,
            bamFile: File
    ): Double {
        val exons = mutableMapOf<Int, Int>()
        val junctions = mutableMapOf<Int, Int>()

        for ( i in (regionStart - 20)..regionStart ) {
            exons[i] = 0
        }

        for ( i in (regionEnd..regionEnd + 20) ) {
            exons[i] = 0
        }

        for ( i in (regionStart + 1)..(regionEnd - 1) ) {
            junctions[i] = 0
        }

        val tmpReader =  SamReaderFactory
                .makeDefault()
                .open(bamFile)
                .query(chromosome, regionStart - 20, regionEnd + 20, false)

        for ( i in tmpReader ) {

            val sites = this.extractSpliceFromCigar(i)
            for ( j in 0..(sites.size - 2) step 2 ) {
                for ( k in sites[j]..sites[j + 1] ) {
                    if ( exons.containsKey(k) ) {
                        exons[k] = 1 + exons[k]!!
                    }

                    if ( junctions.containsKey(k) ) {
                        junctions[k] = 1 + junctions[k]!!
                    }
                }
            }
        }

        fun getMean(data: Map<Int, Int>): Double {
            return data.values.sum().toDouble() / data.size
        }

        return getMean( junctions ) / getMean( exons )
    }
}


