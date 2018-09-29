package com.splice4k.index

import com.splice4k.base.Genes
import com.splice4k.base.SpliceGraph
import com.splice4k.progressbar.ProgressBar
import com.splice4k.tools.FileValidator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.io.PrintWriter
import java.util.*
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since ???
 * @version 20180929
 */


/**
 * 读取提取Splice4k得到的SJ文本格式，或者进一步支持STAR的out.tab格式
 * @param infile 输入文件
 * @param filter 过滤低丰度的SJ
 * @param smrt 是否为SMRT算法运行
 * @param silent 是否print出详细信息
 */
class SJIndex(
        val infile: File,
        val filter: Int,
        val silent: Boolean,
        val smrt: Boolean = false
) {
    val fileFormat: String = FileValidator().check(this.infile)
    private val logger = Logger.getLogger(SJIndex::class.java)
    val transcripts = mutableListOf<Genes>()
    val data = mutableMapOf<String, SpliceGraph>()

    init {
        if ( smrt && this.fileFormat != "bam" ) {
            this.logger.error("Please check input format")
            exitProcess(2)
        }

        this.getAllSJ()
    }


    /**
     * 获取所有的可变剪接事件
     */
    private fun getAllSJ() {

        when ( fileFormat ) {
            "sj" -> this.readSJ( star = false )
            "star" -> this.readSJ( star = true )
            "bam" -> this.readBam()
            else -> {
                this.logger.error("Please check input file format")
                exitProcess(2)
            }
        }

        if ( this.data.isEmpty() ) {
            this.logger.error("${this.infile} format error")
            exitProcess("format error".toInt())
        }
    }


    /**
     * 从extracted splice junctions或者STAR SJ.out.tab文件读取剪接事件
     */
    private fun readSJ( star: Boolean ) {
        logger.info("Reading from ${this.infile}")

        val pattern = when ( star ) {
            false -> "^([\\w\\.]+):(\\d+)-(\\d+)([+-\\.]?)\t(\\d+)$".toRegex()
            else -> "^([\\w\\.]+)\\s(\\d+)\\s(\\d+)\\s([12])\\s\\d+\\s\\d+\\s+(\\d+).*$".toRegex()
        }

        val reader = Scanner(this.infile)

        val pb = ProgressBar(message = "Reading SJ")

        while (reader.hasNext()) {
            val line = reader.nextLine()
            pb.step()
            try{
                var (chromosome, tmpStart, tmpEnd, strand, count) = pattern.find(line)!!.destructured

                strand = when( strand  ) {
                    "" -> "."
                    "1" -> "+"
                    "2" -> "-"
                    else -> strand
                }

                if ( count.toInt() < this.filter ) {
                    continue
                }

                val key = "$chromosome$strand"
                val tmpGraph = when ( this.data.containsKey(key) ) {
                    true -> this.data[key]!!
                    else -> SpliceGraph(chromosome = chromosome, strand = strand.toCharArray().first())
                }

                tmpGraph.addEdge(start = tmpStart.toInt(), end = tmpEnd.toInt(), freq = count.toInt())

                this.data[key] = tmpGraph
            } catch ( e: NullPointerException ) {
                continue
            }
        }
        pb.close()
    }


    /**
     * private function
     * 从单条bam中提取cigar中的N区域
     * @param record 单条SAM/BAM的信息
     * @return 列表，记录了所有的intron的边界信息
     */
    private fun extractSpliceFromCigar( record: SAMRecord): List<Int> {
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

                    results.add(position - 1)
                }
                tmp.clear()
            }
        }

        return results
    }


    /**
     * 从Bam文件中提取所有Splice junctions和其转录本等信息
     */
    private fun readBam() {
        val tmpReader =  SamReaderFactory
                .makeDefault()
                .open(this.infile)
                .iterator()

        this.logger.info("Reading from ${this.infile}")
        val pb = ProgressBar(message = "Reading from Bam")
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

            // SGS构建junctions map
            val tmpGraph = when( this.data.containsKey("${record.referenceName}$strand")  ) {
                true -> this.data["${record.referenceName}$strand"]!!
                else -> SpliceGraph(record.referenceName, strand)
            }

            for ( i in 0..(spliceSites.size - 1) step 2) {
                tmpGraph.addEdge(start = spliceSites[i], end = spliceSites[i + 1])
            }

            this.data["${record.referenceName}$strand"] = tmpGraph

            // SMRT 构建list of transcripts
            if ( smrt ) {
                // 判断reads是否为unique mapped
                if (record.hasAttribute("NH")) {

                    val mapped = record.getAttribute("NH").toString().toInt()

                    if ( mapped <= 0 || mapped > 1) continue
                } else {
                    // 没有NH标签的reads，通常也会造成其他错误，因此直接放弃
                    if (!this.silent) this.logger.warn("${record.readName} does not have attribute NH")
                    continue
                }


                val junctions = this.extractSpliceFromCigar(record)
                // init Genes
                val tmpGene = Genes(
                        chromosome = record.referenceName,
                        start = record.alignmentStart,
                        end = record.alignmentEnd,
                        geneName = record.readName,
                        strand = strand
                )
                // construct exon to transcripts
                tmpGene.exons.addAll(junctions)

                this.transcripts.add(tmpGene)
            }
        }
        pb.close()
    }


    /**
     * 将提取到的Junctions输出到文件
     * @param output 输入文件
     */
    fun writeTo( output: File ) {

        try{
            if (!output.absoluteFile.parentFile.exists()) {
                output.absoluteFile.parentFile.mkdirs()
            }

            val writer = PrintWriter(output)

            for ( v in this.data.values ) {
                writer.print(v.toString())
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