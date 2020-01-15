package com.splice4k.tools

import htsjdk.samtools.*
import java.io.File
import java.nio.BufferUnderflowException

/**
 * 专门用于IR事件中PSI的计算
 * @author Zhang Yiming
 * @since 2018.09.29
 * @version 20180930
 *
 * 20180930 支持在没有index的情况下，自动生成index
 */



class PsiOfIR {

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
            bamFile: File?
    ): Double? {
        if ( bamFile == null ) {
            return null
        }

        val exons = mutableMapOf<Int, Int>()
        val junctions = mutableMapOf<Int, Int>()

        var tmpReader =  SamReaderFactory
                .makeDefault()
                .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                .validationStringency(ValidationStringency.LENIENT)
                .open(bamFile)
        // bai index needed, if .bai not found generate one
        try{
            if ( !File("$bamFile.bai").exists() ) {

                println("Creating index for $bamFile")
                BAMIndexer.createIndex(tmpReader, File("$bamFile.bai"))

                tmpReader.close()
            }
        } catch ( e: htsjdk.samtools.SAMException ) {
            println("Create index failed for $bamFile, ${e.localizedMessage}")
        } finally {
            tmpReader.close()
        }


        // read from BAM/SAM file
        tmpReader =  SamReaderFactory
                .makeDefault()
                .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                .validationStringency(ValidationStringency.LENIENT)
                .open(bamFile)
        try{

            tmpReader.use {
                for ( i in tmpReader.query( chromosome, regionStart - 20, regionEnd + 20, false ) ) {

                    val sites = this.extractSpliceFromCigar(i)
                    for ( j in 0..(sites.size - 2) step 2 ) {

                        val startSite = kotlin.math.max(sites[j], regionStart - 20)
                        val endSite = kotlin.math.min(sites[j + 1], regionEnd + 20)

                        for ( k in startSite..endSite ) {
                            when ( k in regionStart..regionEnd) {
                                true -> junctions[k] = junctions[k] ?: 0 + 1
                                false -> exons[k] = exons[k] ?: 0 + 1
                            }
                        }
                    }
                }
            }

        } catch (e: SAMException) {
            return null
        } catch (e: BufferUnderflowException) {
            println("Can't get reads around $chromosome:$regionStart-$regionEnd from ${bamFile.name}")
            return null
        } finally {
            tmpReader.close()
        }

        fun getMean( data: Map<Int, Int> ): Double {
            return data.values.sum().toDouble() / data.size
        }

        return getMean( junctions ) / getMean( exons )
    }
}


