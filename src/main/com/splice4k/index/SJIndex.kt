package com.splice4k.index

import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.base.JunctionsGraph
import com.splice4k.tools.FileValidator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import me.tongfei.progressbar.ProgressBar
import me.tongfei.progressbar.ProgressBarBuilder
import me.tongfei.progressbar.ProgressBarStyle
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
    private val logger = Logger.getLogger(SJIndex::class.java)
    val transcripts = mutableListOf<Genes>()
    val fileFormat = FileValidator().check( this.infile )
    // chromosome and splice graph
    val data = mutableMapOf<String, JunctionsGraph>()

    init {
        if ( smrt && this.fileFormat !in arrayOf("bam", "gmap") ) {
            this.logger.error("Please check ${this.infile} format, it should be BAM|SAM or extract file from gmap")
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
            "gmap" -> this.readGmap()
            else -> {
                this.logger.error("Please check ${this.infile} format")
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
        val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
        val reader = Scanner(ProgressBar.wrap(FileInputStream(this.infile), pbb))

        reader.use {
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
                        else -> JunctionsGraph(chromosome = info["chromosome"]!!, strand = info["strand"]!!.toCharArray().first())
                    }

                    tmpGraph.addEdge(start = info["start"]!!.toInt(), end = info["end"]!!.toInt(), freq = info["count"]!!.toInt())

                    this.data[key] = tmpGraph
                } catch ( e: NullPointerException ) {
                    continue
                }
            }
        }

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

                // 'D', 'H', 'P',
                if (i !in arrayOf('S', 'D') ) {
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
        val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
        val pb = ProgressBar.wrap(tmpReader, pbb)
        val junctions = mutableMapOf<String, Int>()

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


            // 统计所有的junctions的数量，便于filter
            for ( i in 1..(spliceSites.size - 2) step 2) {
                val key = "${record.referenceName}\t${spliceSites[i]}\t${spliceSites[i + 1]}\t$strand"
                junctions[key] = junctions[key]?: 0 + 1
            }

            // SMRT 构建list of transcripts
            if ( this.smrt ) {

                // init Genes
                val tmpGene = Genes(
                        chromosome = record.referenceName,
                        start = record.alignmentStart,
                        end = record.alignmentEnd,
                        geneName = record.readName,
                        strand = strand
                )

                // construct exon to transcripts
                try{
                    for ( i in 0..(spliceSites.size - 2) step 2 ) {
                        tmpGene.exons.add(Exons(
                                chromosome = record.referenceName,
                                start = spliceSites[i],
                                end = spliceSites[i + 1],
                                exonId = ""
                        ))
                    }
                } catch ( e: com.splice4k.errors.ChromosomeException ) {
                    this.logger.error(e.localizedMessage)
                    println(record.readName)
                    println("${record.start}\t${record.end}")
                    println(spliceSites)
                    exitProcess(0)
                }

                this.transcripts.add(tmpGene)
            }
        }

        tmpReader.close()

        for ( (key, v) in junctions ) {
            if ( v < this.filter ) {
                continue
            }

            val intron = key.split("\t")

            val dataKey = "${intron[0]}${intron[3]}"

            this.data[dataKey] = this.data[dataKey] ?: JunctionsGraph(intron[0], intron[3].toCharArray()[0])

            this.data[dataKey]!!.addEdge(intron[1].toInt(), intron[2].toInt(), freq = v)
        }
    }


    /**
     * 从gmap align的文件中提取读取reads来做三代
     */
    private fun readGmap() {

        this.logger.info("Reading from ${this.infile}")
        val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
        val reader = Scanner(ProgressBar.wrap(FileInputStream(this.infile), pbb))
        val junctions = mutableMapOf<String, Int>()


        while ( reader.hasNext() ) {
            val line = reader.nextLine()
            val lines = line.split("\t")

            if ( lines[4] == "" ) {
                continue
            }

            val chromosome = lines[0]
            val start = lines[1].toInt()
            val end = lines[2].toInt()
            val strand = lines[3].toCharArray()[0]

            val introns = lines[4].split(",").map { it.toInt() }
            val exons = mutableListOf(start)

            for ( i in 0..(introns.size - 2) step 2) {
                val key = "$chromosome\t${introns[i]}\t${introns[i + 1]}\t$strand"

                junctions[key] = junctions[key]?: 0 + 1

                // construct exons
                exons.add(introns[i] - 1)
                exons.add(introns[i + 1] + 1)
            }

            exons.add(end)

            // init Genes
            val tmpGene = Genes(
                    chromosome = chromosome,
                    start = start,
                    end = end,
                    geneName = "NA",
                    strand = strand
            )

            // construct exon to transcripts
            try{
                for ( i in 0..(exons.size - 2) step 2 ) {
                    tmpGene.exons.add(Exons(
                            chromosome = chromosome,
                            start = exons[i],
                            end = exons[i + 1],
                            exonId = ""
                    ))
                }
            } catch ( e: com.splice4k.errors.ChromosomeException ) {
                this.logger.error(e.localizedMessage)
                println(exons)
                exitProcess(0)
            }

            this.transcripts.add(tmpGene)
        }

        reader.close()

        for ( (key, v) in junctions ) {
            if ( v < this.filter ) {
                continue
            }

            val intron = key.split("\t")

            val dataKey = "${intron[0]}${intron[3]}"

            this.data[dataKey] = this.data[dataKey] ?: JunctionsGraph(intron[0], intron[3].toCharArray()[0])

            this.data[dataKey]!!.addEdge(intron[1].toInt(), intron[2].toInt(), freq = v)
        }

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