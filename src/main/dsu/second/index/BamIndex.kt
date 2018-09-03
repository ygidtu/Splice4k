package dsu.second.index

import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMRecord
import java.io.File
import org.apache.log4j.Logger

import dsu.carrier.GenomicLoci
import dsu.progressbar.ProgressBar



/**
  * @author Zhangyiming
  * @since 2018.09.03
  * @version 20180903
  * Extract Splice Junctions from Bam/Sam file
  */


class BamIndex(bam: String) {
    private val logger: Logger = Logger.getLogger(BamIndex::class.java)
    private val reader = SamReaderFactory
            .makeDefault()
            .open(File(bam))

    val data: Map<GenomicLoci, Int>

    init {
        data = this.getAllSJ()
    }

    /**
     * private function
     * 从单条bam中提取cigar中的N区域
     * @param record 单条SAM/BAM的信息
     * @return 列表，记录了所有的intron的边界信息
     */
    private fun extractSpliceFromCigar(record: SAMRecord): List<Int> {
        val results: MutableList<Int> = mutableListOf()
        var position = record.alignmentStart
        val tmp = mutableListOf<Char>()

        for (i in record.cigar.toString()) {
            if (i in '0'..'9') {  // 如果是数字，就加到list中
                tmp.add(i)
            } else {
                if (tmp.size == 0) {
                    continue
                }

                // Soft clip以及insertion的两种区域都不记载在alignment之内
                if (i != 'S' && i != 'I') {
                    position += tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                }

                if (i == 'N') {
                    results.add(
                            position - tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                    )

                    results.add(position)
                }
                tmp.clear()
            }
        }

        return results
    }

    /**
     * 获取bam文件中的所有SJ及其freq数值
     * @return Map<GenomicLoci, Int>
     */
    private fun getAllSJ(): Map<GenomicLoci, Int> {
        val tmpReader = this.reader.iterator()

        val results = mutableMapOf<GenomicLoci, Int>()

        logger.info("Star Reading Splice Junctions from Bam")
        val pb = ProgressBar(message = "Reading SJ from Bam")
        for ( record in tmpReader) {

            pb.step()

            val spliceSites = this.extractSpliceFromCigar(record)

            for ( i in 0..spliceSites.size step 2) {
                val tmpLoci = GenomicLoci(
                        chromosome = record.referenceName,
                        start = spliceSites[i],
                        end = spliceSites[i + 1]
                )

                when {
                    results.containsKey(tmpLoci) -> results[tmpLoci] = results[tmpLoci]!! + 1
                    else -> results[tmpLoci] = 0
                }
            }
        }

        return results
    }

 }