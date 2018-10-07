package com.splice4k.index

import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.base.SpliceGraph
import com.splice4k.tools.FileValidator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import me.tongfei.progressbar.ProgressBar
import org.apache.log4j.Logger
import java.io.File
import java.io.FileInputStream
import java.io.IOException
import java.io.PrintWriter
import java.util.*
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since ???
 * @version 20181006
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
        private val filter: Int,
        private val silent: Boolean,
        private val smrt: Boolean = false
) {
    val fileFormat: String = FileValidator().check(this.infile)
    private val logger = Logger.getLogger(SJIndex::class.java)
    val transcripts = mutableListOf<Genes>()

    // chromosome and splice graph
    val data = mutableMapOf<String, SpliceGraph>()

    init {
        if ( smrt && this.fileFormat != "bam" ) {
            this.logger.error("Please check ${this.infile} format, it should be BAM|SAM")
            exitProcess(2)
        }

        this.getAllSJ()
    }


    /**
     * 获取所有的可变剪接事件
     */
    private fun getAllSJ() {

        when ( fileFormat ) {
            "sj" -> this.readSJ()
            "star" -> this.readSJ()
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


    /*
        val pattern = when ( star ) {
            false -> "^([\\w\\.]+):(\\d+)-(\\d+)([+-\\.]?)\t(\\d+)$".toRegex()
            else -> when (this.unique) {
            true -> "^([\\w\\.]+)\\s(\\d+)\\s(\\d+)\\s([12])\\s\\d+\\s1\\s+(\\d+).*$".toRegex()
            else -> "^([\\w\\.]+)\\s(\\d+)\\s(\\d+)\\s([12])\\s\\d+\\s[01]\\s+(\\d+).*$".toRegex()
            }
        }

        var (chromosome, tmpStart, tmpEnd, strand, count) = pattern.find(line)!!.destructured

        strand = when( strand  ) {
            "" -> "."
            "1" -> "+"
            "2" -> "-"
            else -> strand
        }
     */


    /**
     * 从sj文件中通过split的方式获取数据。
     * @param line 文件行
     * @return 关键字map
     */
    private fun getSitesFromSJ(line: String): Map<String, String> {
        val results = mutableMapOf<String, String>()

        var lines = line.replace("\n", "").split("\t")

        results["count"] = lines.last()

        results["strand"] = when {
            lines.first().endsWith("+") -> "+"
            lines.first().endsWith("-") -> "-"
            else -> "."
        }

        lines = lines.first().split(":")

        results["chromosome"] = lines.first()

        lines = lines.last().split("-")

        results["start"] = lines.first()

        results["end"] = lines[1].replace("[+-\\.]".toRegex(), "")

        return results
    }


    /**
     * 从star文件中通过split的方式获取数据。
     * @param line 文件行
     * @return 关键字map
     */
    private fun getSitesFromSTAR(line: String): Map<String, String> {
        val lines = line.split("\\s+".toRegex())

        return mapOf(
                "chromosome" to lines[0],
                "start" to lines[1],
                "end" to lines[2],
                "strand" to when (lines[3]) {
                    "1" -> "+"
                    "2" -> "-"
                    else -> "."
                },
                "uniq" to lines[5],
                "count" to lines[6]
        )
    }


    /**
     * 从extracted splice junctions或者STAR SJ.out.tab文件读取剪接事件
     */
    private fun readSJ() {
        logger.info("Reading from ${this.infile}")

        val reader = Scanner(ProgressBar.wrap(FileInputStream(this.infile), "Reading"))

        while (reader.hasNext()) {
            val line = reader.nextLine()

            try{

                val info = when(this.fileFormat) {
                    "sj" -> getSitesFromSJ(line)
                    else -> getSitesFromSTAR(line)
                }

                if ( info["count"]!!.toInt() < this.filter ) {
                    continue
                }

                val key = "${info["chromosome"]!!}${info["strand"]!!}"
                val tmpGraph = when ( this.data.containsKey(key) ) {
                    true -> this.data[key]!!
                    else -> SpliceGraph(chromosome = info["chromosome"]!!, strand = info["strand"]!!.toCharArray().first())
                }

                tmpGraph.addEdge(start = info["start"]!!.toInt(), end = info["end"]!!.toInt(), freq = info["count"]!!.toInt())

                this.data[key] = tmpGraph
            } catch ( e: NullPointerException ) {
                continue
            }
        }
        reader.close()
    }


    /**
     * private function
     * 从单条bam中提取cigar中的N区域
     * @param record 单条SAM/BAM的信息
     * @return 列表，记录了所有的intron的边界信息
     */
    private fun extractSpliceFromCigar( record: SAMRecord ): List<Int> {
        var position = record.alignmentStart
        val tmp = mutableListOf<Char>()
        val results: MutableList<Int> = mutableListOf( position )

        for ( i in record.cigarString ) {
            if (i in '0'..'9') {  // 如果是数字，就加到list中
                tmp.add(i)
            } else {
                if (tmp.size == 0) {
                    continue
                }

                if (i != 'D' && i != 'H' && i != 'P' ) {
                    position += tmp.joinToString(separator = "").toInt()
                }

                if (i == 'N') {
                    results.add(
                            position - tmp.joinToString(separator = "").toInt()
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
     * 从Bam文件中提取所有Splice junctions和其转录本等信息
     */
    private fun readBam() {
        val tmpReader =  SamReaderFactory
                .makeDefault()
                .open(this.infile)
                .iterator()

        this.logger.info("Reading from ${this.infile}")

        val pb = ProgressBar.wrap(tmpReader.stream(), "Reading")

        for ( record in pb) {

            // 判断reads是否为unique mapped
            if (record.hasAttribute("NH")) {

                val mapped = record.getAttribute("NH").toString().toInt()

                if ( mapped <= 0 || mapped > 1) continue
            } else {
                // 没有NH标签的reads，通常也会造成其他错误，因此直接放弃
                if (!this.silent) this.logger.warn("${record.readName} does not have attribute NH")
                continue
            }


            val spliceSites = this.extractSpliceFromCigar(record)


            /*
             这个问题以前从没注意过，跟STAR比对过才发现在这个问题
             mate reverse strand或者read reverse strand都有可能代表了这条reads真实存在的链情况
             STAR的解决方案可能是按照first in pair来处理的，
             只要是first in pair，就看是mate的情况；否则看reads negate strand
             其实都代表这条reads在负链上
              */
            val strandJudge = when( record.readPairedFlag ) {
                true -> when( record.firstOfPairFlag ) {
                    true -> record.mateNegativeStrandFlag
                    false -> record.readNegativeStrandFlag
                }

                false -> record.readNegativeStrandFlag
            }

            val strand = when ( strandJudge ) {
                true -> '-'
                false -> '+'
            }


            // SGS构建junctions map
            val tmpGraph = when( this.data.containsKey("${record.referenceName}$strand")  ) {
                true -> this.data["${record.referenceName}$strand"]!!
                else -> SpliceGraph(record.referenceName, strand)
            }

            for ( i in 1..(spliceSites.size - 2) step 2) {
                tmpGraph.addEdge(start = spliceSites[i], end = spliceSites[i + 1])
            }


            this.data["${record.referenceName}$strand"] = tmpGraph

            // SMRT 构建list of transcripts
            if ( smrt ) {

                // init Genes
                val tmpGene = Genes(
                        chromosome = record.referenceName,
                        start = record.alignmentStart,
                        end = record.alignmentEnd,
                        geneName = record.readName,
                        strand = strand
                )

                // construct exon to transcripts
                for ( i in 0..(spliceSites.size - 2) step 2 ) {
                    tmpGene.exons.add(Exons(
                            chromosome = record.referenceName,
                            start = spliceSites[i],
                            end = spliceSites[i + 1],
                            exonId = ""
                    ))
                }

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