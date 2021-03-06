package com.splice4k.index


import com.splice4k.base.Genes
import com.splice4k.base.Reads
import com.splice4k.base.JunctionsGraph
import com.splice4k.tools.FileValidator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import me.tongfei.progressbar.ProgressBar
import me.tongfei.progressbar.ProgressBarBuilder
import me.tongfei.progressbar.ProgressBarStyle
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
 *
 * 2018.10.18 添加对gmap align文件的支持，重新format一下代码结构
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

    val transcripts = mutableListOf<Genes>()
    val reads = mutableListOf<Reads>()
    val fileFormat = FileValidator().check( this.infile )
    // chromosome and splice graph
    val data = mutableMapOf<String, Int>()

    init {
        if ( this.smrt && this.fileFormat !in arrayOf("bam", "gmap", "gmapE") ) {
            println("Please check ${this.infile} format, it should be BAM|SAM or gmap align file")
            exitProcess(2)
        }
        this.getAllSJ()
    }


    /**
     * 获取所有的可变剪接事件
     */
    private fun getAllSJ() {

        println("Reading from ${this.infile}")

        when ( this.fileFormat ) {
            "bam" -> this.readBam( reads = this.reads )
            in arrayOf("gmap", "gmapE") -> this.readGmap( reads = this.reads )
            in arrayOf( "star", "sj" ) -> this.readSJ()
            else -> {
                println("Please check ${this.infile} format")
                exitProcess(2)
            }
        }

        // 构建转录本
        if ( this.smrt ) {
            this.reads.forEach {
                // init Genes
                val tmpGene = Genes(
                        chromosome = it.chromosome,
                        start = it.start,
                        end = it.end,
                        geneName = "NA",
                        strand = it.strand
                )

                // construct exon to transcripts
                try{
                    tmpGene.exons = it.getExon().toMutableList()
                } catch ( e: com.splice4k.errors.ChromosomeException ) {
                    println(e.localizedMessage)
                    println(it.getExon())
                    exitProcess(0)
                }

                this.transcripts.add(tmpGene)
            }
        }
    }


    /**
     * 从extracted splice junctions或者STAR SJ.out.tab文件读取剪接事件
     */
    private fun readSJ() {

        val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
        val reader = Scanner(ProgressBar.wrap(FileInputStream(this.infile), pbb))

        reader.use {
            while (reader.hasNext()) {
                val line = reader.nextLine()

                try{

                    when(this.fileFormat) {
                        "sj" -> {
                            val cleanedLine = line.replace("[:-]".toRegex(), "\t").split("\t")
                            this.data[cleanedLine.subList(0, cleanedLine.size - 1).joinToString("\t")] = cleanedLine.last().toInt()
                        }
                        else -> {
                            val lines = line.split("\\s+".toRegex())
                            val strand = when (lines[3]) {
                                "1" -> "+"
                                "2" -> "-"
                                else -> "."
                            }
                            this.data["${lines[0]}\t${lines[1]}\t${lines[2]}\t$strand"] = lines[6].toInt()
                        }
                    }


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
        val results: MutableList<Int> = mutableListOf()

        for ( i in record.cigarString ) {
            if (i in '0'..'9') {  // 如果是数字，就加到list中
                tmp.add(i)
            } else {
                if (tmp.size == 0) {
                    continue
                }

                if (i in arrayOf('M', 'I', 'N') ) {
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

        results.add(position - 1)
        return results
    }


    /**
     * 从Bam文件中提取所有Splice junctions和其转录本等信息
     * @param junctions junctions的count map
     * @param reads Reads的列表
     */
    private fun readBam( reads: MutableList<Reads> ) {
        val tmpReader =  SamReaderFactory
                .makeDefault()
                .open(this.infile)
                .iterator()

        val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
        val pb = ProgressBar.wrap(tmpReader, pbb)


        for ( record in pb) {

            // 判断reads是否为unique mapped
            if (record.hasAttribute("NH")) {

                val mapped = record.getAttribute("NH").toString().toInt()

                if ( mapped <= 0 || mapped > 1) continue
            } else {
                // 没有NH标签的reads，通常也会造成其他错误，因此直接放弃
                if (!this.silent) println("${record.readName} does not have attribute NH")
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

            reads.add(Reads(
                    chromosome = record.referenceName,
                    start = record.alignmentStart,
                    end = spliceSites.last(),
                    name = record.readName,
                    strand = strand,
                    introns = spliceSites.slice(0..(spliceSites.size - 2))
            ))
        }

        tmpReader.close()
    }


    /**
     * 从gmap align的文件中提取读取reads来做三代
     * @param junctions junctions的count map
     * @param reads Reads的列表
     */
    private fun readGmap( reads: MutableList<Reads> ) {

        val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
        val reader = Scanner(ProgressBar.wrap(FileInputStream(this.infile), pbb))
        val tmpSites = mutableListOf<Int>()     // gmap align extracted sites from single reads
        var tmpStrand: Char? = null                            // gmap align extracted strand from single reads
        var tmpChromo: String? = null                          // same as tmpStrand

        while ( reader.hasNext() ) {
            val line = reader.nextLine()

            if ( this.fileFormat == "gmapE" ) {
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

                    this.data[key] = (this.data[key]?: 0) + 1

                    // construct exons
                    exons.add(introns[i] - 1)
                    exons.add(introns[i + 1] + 1)
                }

                exons.add(end)

                reads.add(Reads(
                        chromosome = chromosome,
                        start = start,
                        end = end,
                        name = "",
                        strand = strand,
                        introns = introns
                ))
            } else if ( this.fileFormat == "gmap" ) {
                val gmapPattern = "^\\s+([+-])([\\w\\.]+):(\\d+)-(\\d+)\\s+\\(\\d+-\\d+\\)\\s+\\d{0,3}%.*"

                if ( line.matches(gmapPattern.toRegex()) ) {
                    val (strand, chromosome, start, end, _) = gmapPattern.toRegex().matchEntire(line)!!.destructured
                    if ( tmpStrand == null ) {
                        tmpStrand = strand.toCharArray()[0]
                    }

                    if ( tmpChromo == null ) {
                        tmpChromo = chromosome
                    }

                    tmpSites.add( start.toInt() )
                    tmpSites.add( end.toInt() )
                } else {
                    if ( tmpSites.size > 2 ) {
                        tmpSites.sort()
                        val introns = mutableListOf<Int>()

                        for ( i in 1..(tmpSites.size - 2) step 2 ) {
                            introns.add( tmpSites[i] + 1 )
                            introns.add( tmpSites[i + 1] - 1 )

                            val key = "$tmpChromo\t${tmpSites[i] + 1}\t${tmpSites[i + 1] - 1}\t$tmpStrand"

                            this.data[key] = (this.data[key]?: 0) + 1
                        }

                        reads.add(Reads(
                                chromosome = tmpChromo!!,
                                start = tmpSites.first(),
                                end = tmpSites.last(),
                                name = tmpChromo,
                                strand = tmpStrand!!,
                                introns = introns
                        ))
                    }

                    tmpChromo = null; tmpStrand = null; tmpSites.clear()
                }
            }

        }

        reader.close()
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
            println(err.message)
            for (i in err.stackTrace) {
                println(i)
            }
        }
    }


    /**
     * get junctionGraph
     * @return List<JunctionGraph>
     */
    fun getJunctionGraph(filtered: HashSet<String>?): List<JunctionsGraph> {
        return this.getJunctionGraphMap(filtered).values.toList()
    }


    /**
     * get junctionGraph separated by chromosome and strands
     * @return Map<String, JunctionGraph>, keys are chromosome+strand
     */
    fun getJunctionGraphMap(filtered: HashSet<String>?): Map<String, JunctionsGraph> {
        val res = mutableMapOf<String, JunctionsGraph>()
        // 构建junction图
        for ( (key, v) in this.data ) {
            if ( v < this.filter ) {
                continue
            }

            if ( filtered != null && key in filtered ) {
                continue
            }

            val intron = key.split("\t")

            val dataKey = "${intron[0]}${intron[3]}"

            res[dataKey] = res[dataKey] ?: JunctionsGraph(intron[0], intron[3].toCharArray()[0])

            res[dataKey]!!.addEdge(intron[1].toInt(), intron[2].toInt(), freq = v)
        }
        return res
    }
}