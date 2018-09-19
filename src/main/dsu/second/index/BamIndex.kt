package dsu.second.index

import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMRecord
import java.io.File

import dsu.progressbar.ProgressBar
import dsu.carrier.SpliceGraph
import java.io.IOException
import java.io.PrintWriter


/**
  * @author Zhangyiming
  * @since 2018.09.03
  * @version 20180903
  * Extract Splice Junctions from Bam/Sam file
  */


class BamIndex(bam: String): SJIndex(bam) {

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
    override fun getAllSJ() {
        val tmpReader =  SamReaderFactory
                .makeDefault()
                .open(File(this.infile))
                .iterator()

        logger.info("Start Reading Splice Junctions from Bam")
        val pb = ProgressBar(message = "Reading SJ from Bam")
        for ( record in tmpReader) {

            pb.step()

            val spliceSites = this.extractSpliceFromCigar(record)

            if (spliceSites.isEmpty()) {
                continue
            }

            val strand = when(record.readNegativeStrandFlag) {
                true -> '-'
                false -> '+'
            }

            val tmpGraph = when( this.data.containsKey("${record.referenceName}$strand")  ) {
                true -> this.data["${record.referenceName}$strand"]!!
                else -> SpliceGraph(record.referenceName, strand)
            }

            for ( i in 0..(spliceSites.size - 1) step 2) {
                tmpGraph.addEdge(start = spliceSites[i] - 1, end = spliceSites[i + 1] + 1)
            }

            this.data["${record.referenceName}$strand"] = tmpGraph
        }
    }


    /**
     * 将提取到的Junctions输出到文件
     * @param output 输入文件
     */
    fun writeTo(output: File) {

        try{
            if (!output.absoluteFile.parentFile.exists()) output.absoluteFile.parentFile.mkdirs()

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