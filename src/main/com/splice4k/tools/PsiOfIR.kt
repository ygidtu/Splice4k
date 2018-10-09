package com.splice4k.tools

import htsjdk.samtools.*
import org.apache.log4j.Logger
import java.io.File

/**
 * 专门用于IR事件中PSI的计算
 * @author Zhang Yiming
 * @since 2018.09.29
 * @version 20180930
 *
 * 20180930 支持在没有index的情况下，自动生成index
 */



class PsiOfIR {
    private val logger = Logger.getLogger(PsiOfIR::class.java)

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
    ): Double? {
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


        try{
            val indexFile = File("${bamFile.absolutePath}.bai")
            if ( !indexFile.exists() ) {
                val tmpReader =  SamReaderFactory
                        .makeDefault()
                        .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                        .validationStringency(ValidationStringency.LENIENT)
                        .open(bamFile)

                this.logger.info("Creating index for $bamFile")
                BAMIndexer.createIndex(tmpReader, indexFile)

                tmpReader.close()
            }

            val tmpReader =  SamReaderFactory
                    .makeDefault()
                    .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                    .validationStringency(ValidationStringency.LENIENT)
                    .open(bamFile)


            for ( i in tmpReader.query( chromosome, regionStart, regionEnd, false ) ) {

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

            tmpReader.close()
        } catch (e: SAMException) {
            return null
        }

        fun getMean(data: Map<Int, Int>): Double {
            return data.values.sum().toDouble() / data.size
        }

        return getMean( junctions ) / getMean( exons )
    }
}


