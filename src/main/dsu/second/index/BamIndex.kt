package dsu.second.index

import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMRecord
import java.io.File
import org.apache.log4j.Logger

import dsu.carrier.Exons
import dsu.progressbar.ProgressBar
import java.io.IOException
import java.io.PrintWriter
import java.util.*


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

    val data = mutableMapOf<Exons, Int>()
    val sameStart = mutableMapOf<Int, MutableList<Exons>>()
    val sameEnd = mutableMapOf<Int, MutableList<Exons>>()

    init {
        this.getAllSJ()
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

    private fun getSite(loci: Exons, start: Boolean = true): Int {
        return when(start) {
            true -> Objects.hash(loci.chromosome, loci.start, loci.strand)
            false -> Objects.hash(loci.chromosome, loci.end, loci.strand)
        }
    }

    /**
     * 获取bam文件中的所有SJ及其freq数值
     * @return Map<GenomicLoci, Int>
     */
    private fun getAllSJ() {
        val tmpReader = this.reader.iterator()

        logger.info("Star Reading Splice Junctions from Bam")
        val pb = ProgressBar(message = "Reading SJ from Bam")
        for ( record in tmpReader) {

            pb.step()

            val spliceSites = this.extractSpliceFromCigar(record)

            for ( i in 0..spliceSites.size step 2) {
                val tmpLoci = Exons(
                        chromosome = record.referenceName,
                        start = spliceSites[i],
                        end = spliceSites[i + 1],
                        strand = when (record.readNegativeStrandFlag) {
                            true -> '-'
                            else -> '+'
                        }
                )

                when {
                    this.data.containsKey(tmpLoci) -> this.data[tmpLoci] = this.data[tmpLoci]!! + 1
                    else -> this.data[tmpLoci] = 0
                }

                val start = getSite(tmpLoci, true)
                val end = getSite(tmpLoci, false)
                when (this.sameStart.containsKey(start)) {
                    true -> this.sameStart[start]!!.add(tmpLoci)
                    false -> this.sameStart[start] = mutableListOf(tmpLoci)
                }

                when (this.sameEnd.containsKey(end)) {
                    true -> this.sameEnd[end]!!.add(tmpLoci)
                    false -> this.sameEnd[end] = mutableListOf(tmpLoci)
                }
            }
        }
    }


    /**
     * 将提取到的Junctions输出到文件
     * @param output 输入文件
     */
    fun writeTo(output: File) {

        try{
            if (!output.parentFile.exists()) output.parentFile.mkdirs()

            val writer = PrintWriter(output)

            for ( (k, v) in this.data.iterator() ) {
                writer.println("$k\t$v")
            }
            writer.close()
        } catch (err: IOException) {
            logger.error(err.message)
            for (i in err.stackTrace) {
                logger.error(i)
            }
        }
    }

 }