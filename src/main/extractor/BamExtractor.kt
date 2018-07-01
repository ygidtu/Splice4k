package main.extractor

import java.io.File

import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMRecord

import org.apache.log4j.Logger

import main.carrier.Genes
import main.errors.ExonException


/**
 * @since 2018.06.18
 * @version 0.1
 * @author zhangyiming
 *
 * 从bam/sam文件中提取所有的exon区域，并且统计其exon数目
 */


/**
 * @param bam 输入的bam文件的路径
 * @param unique 保留匹配到几个位点上的reads， 默认为1
 * @param silent Boolean值，减少信息输出，默认为false
 */
class BamExtractor(
        bam: String,
        private val unique: Int=1,
        private val silent: Boolean = false
): Extractor(silent) {

    private val logger: Logger = Logger.getLogger(BamExtractor::class.java)
    private val reader = SamReaderFactory
            .makeDefault()
            .open(File(bam))

    init {
        this.data = getAllRegions()
        this.totalLine = this.data.size
    }


    /**
     * private function
     * 从单条bam中提取cigar中的N区域
     * @param record 单条SAM/BAM的信息
     * @return 列表，记录了所有的intron的边界信息
     */
    private fun extractSpliceFromCigar(record: SAMRecord): List<Int> {
        val results = mutableListOf(record.alignmentStart)
        var position = record.alignmentStart
        val tmp = mutableListOf<Char>()

        for (i in record.cigar.toString()) {
            if (i in '0'..'9') {  // 如果是数字，就加到list中
                tmp.add(i)
            } else {
                if (tmp.size == 0) {
                    continue
                }
                position += tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                if (i == 'N') {
                    results.add(
                            position - tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                    )

                    results.add(position)
                }
                tmp.clear()
            }
        }
        results.add(record.alignmentEnd)
        return results
    }


    /**
     * 获取每条bam中的exon信息
     * @return list of Genes
     */
    private fun getAllRegions(): List<Genes> {
        val results = mutableListOf<Genes>()
        var gap = 1

        for ( (index, record) in this.reader.withIndex()) {

            if (index % gap == 0) {
                this.logger.info("Reading $index lines")
                if (gap < 10001) gap *= 10
            }

            // 判断reads是否为unique mapped
            if (record.hasAttribute("NH")) {

                val mapped = record.getAttribute("NH").toString().toInt()

                if ( mapped <= 0 || mapped > this.unique) continue
            } else {
                // 没有NH标签的reads，通常也会造成其他错误，因此直接放弃
                if (!this.silent) this.logger.warn("${record.readName} does not have attribute NH")
                continue
            }

            // get all exons
            val introns = this.extractSpliceFromCigar(record)

            if (introns.size < 2) continue

            // init Genes
            val tmpGene = Genes(
                    chrom = record.referenceName,
                    start = record.alignmentStart,
                    end = record.alignmentEnd,
                    geneName = record.readName,
                    strand = when(record.readNegativeStrandFlag) {
                        true -> '-'
                        false -> '+'
                    }
            )

            // add exons
            val tmpExons = mutableListOf<Array<Int>>()
            for (i in 0..(introns.size - 1) step 2) {
                if (introns[i] > introns[i + 1]) throw ExonException("start[${introns[i]}] > end[${introns[i+1]}]")
                tmpExons.add(arrayOf(introns[i], introns[i + 1]))
            }
            tmpGene.exons = tmpExons
            results.add(tmpGene)
        }
        results.sortWith(compareBy({it.chrom}, {it.start}, {it.end}))
        return results
    }


    fun get(): List<Genes> {
        return this.data
    }
}

/*
fun main(args: Array<String>) {
    /*
    GT-AG规则
    */
    val test = BamExtractor("/home/zhang/splicehunter_test/test.bam")

    test.saveTo("/home/zhang/splicehunter_test/bam_extracted.txt")
}
*/
