package com.splice4k.index

import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.progressbar.ProgressBar
import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.lang.NullPointerException
import java.util.*
import kotlin.system.exitProcess
import com.splice4k.tools.FileValidator


/**
 * 注释文件index的基础类
 * @version 2018.9.29
 * @author Zhang Yiming
 * @since ???
 */


/**
 * @param infile 输入文件
 * @param smrt 是否为SMRT算法运行
 */
class AnnotationIndex(
        private val infile: File,
        val smrt: Boolean = false
) {
    private val logger = Logger.getLogger(AnnotationIndex::class.java)
    private val fileFormat = FileValidator().check(this.infile)
    val data = mutableMapOf<String, MutableList<Exons>>()

    val transcripts: MutableList<Genes> = mutableListOf()

    init {

        when( this.fileFormat ) {
            "gtf" ->  this.readExonFromGtf()

            "gff" -> this.readExonFromGff()

            else -> {
                logger.info("Please check reference file format")
                exitProcess(2)
            }
        }
    }


    /**
     * 从Gff文件中读取信息
     */
    private fun readExonFromGff() {
        fun getSource(info: List<String>): Map<String, String> {
            val results: MutableMap<String, String> = mutableMapOf()

            for (inf in info[0].split(";")) {
                val tmp = inf.split("=")

                if ( ":" in tmp[1] ) {
                    results[tmp[0]] = tmp[1].split(":")[1]
                } else {
                    results[tmp[0]] = tmp[1]
                }
            }

            return results
        }

        try {
            val reader = Scanner(this.infile)

            this.logger.info("Reading from ${this.infile}")
            val pb = ProgressBar(message = "Reading exon from Gff")

            val geneTranscript = mutableMapOf<String, String>()

            val exons: MutableMap<String, MutableList<Exons>> = mutableMapOf()


            while (reader.hasNext()) {
                val line = reader.nextLine()

                if ( line.startsWith("#") ) {
                    continue
                }

                val lines = line.replace("\n", "").split("\t".toRegex())
                pb.step()

                val sources = getSource(lines.subList(8, lines.size))

                if (  // 两种标准，严防不标准的gff文件
                        "transcript_id" in sources.keys ||
                        (lines[2].matches(".*(rna|transcript).*".toRegex(RegexOption.IGNORE_CASE)) &&
                                !lines[2].matches(".*gene.*".toRegex(RegexOption.IGNORE_CASE)))
                ) {

                    geneTranscript[sources["ID"]!!] = sources["Parent"]!!

                    if ( this.smrt ) {
                        this.transcripts.add(Genes(
                                chromosome = lines[0],
                                start = lines[3].toInt(),
                                end = lines[4].toInt(),
                                strand = lines[6].toCharArray()[0],
                                information = sources
                        ))
                    }

                } else if (lines[2] == "exon") {

                    val tmp = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            exonId = sources["exon_id"] ?: sources["ID"]!!
                    )

                    try{

                        tmp.source["transcript"]!!.add(sources["Parent"]!!)
                        tmp.source["gene"]!!.add(geneTranscript[sources["Parent"]]!!)

                    } catch (e: NullPointerException) {
                        println(line)
                        println(sources)
                        exitProcess(0)
                    }


                    val tmpExons = mutableListOf(tmp)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmpExons)
                    } else {
                        this.data[key] = tmpExons
                    }

                    if ( this.smrt ) {
                        val tmpExon = mutableListOf(tmp)

                        if (exons.containsKey(sources["Parent"]!!)) {
                            tmpExon.addAll(exons[sources["Parent"]!!]!!)
                        }
                        exons[sources["Parent"]!!] = tmpExon
                    }
                }
            }

            reader.close()
            pb.close()

            if ( this.smrt ) {
                // 为转录本指定外显子
                for (i in this.transcripts) {
                    if (exons.containsKey(i.transcriptId)) {
                        i.exons.addAll(exons[i.transcriptId]!!)
                    }
                }
            }

        }catch (e: IOException) {
            this.logger.error(e.toString())
        }
    }


    /**
     * 从Gtf文件中读取信息
     */
    private fun readExonFromGtf() {
        fun getSource(info: List<String>): Map<String, String> {
            val results: MutableMap<String, String> = mutableMapOf()

            for (i in 0..(info.size - 1) step 2) {
                try{
                    results[info[i]] = info[i+1]     // get rid of useless characters
                            .replace("\"", "")
                            .replace(";", "")
                } catch ( e: IndexOutOfBoundsException ) {
                    continue
                }


            }
            return results
        }


        try {
            val reader = Scanner(this.infile)

            this.logger.info("Reading from ${this.infile}")
            val pb = ProgressBar(message = "Reading exon from Gtf")

            val exons: MutableMap<String, MutableList<Exons>> = mutableMapOf()
            var exonNumber = 1
            while (reader.hasNext()) {
                val line = reader.nextLine()
                val lines = line.split("\\s+".toRegex())
                pb.step()

                if ( line.startsWith("#") ) {
                    continue
                }

                val sources = getSource(lines.subList(8, lines.size))

                if ( this.smrt && lines[2] == "transcript" ) {
                    this.transcripts.add(Genes(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            information = sources
                    ))

                    exonNumber = 1
                } else if (lines[2] == "exon") {

                    var exonId = sources["exon_id"] ?: sources["ID"]

                    if ( exonId == null ) {
                        exonId = "${sources["transcript_id"]!!}$exonNumber"
                        exonNumber++
                    }

                    val tmp = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            exonId = exonId
                    )

                    tmp.source["transcript"]!!.add(sources["transcript_id"]!!)
                    tmp.source["gene"]!!.add(sources["gene_id"]!!)

                    val tmpExons = mutableListOf(tmp)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmpExons)
                    } else {
                        this.data[key] = tmpExons
                    }

                    if ( this.smrt ) {
                        val tmpExon = mutableListOf(tmp)

                        if (exons.containsKey(sources["transcript_id"]!!)) {
                            tmpExon.addAll(exons[sources["transcript_id"]!!]!!)
                        }
                        exons[sources["transcript_id"]!!] = tmpExon
                    }
                }
            }

            reader.close()
            pb.close()


            if ( this.smrt ) {
                for (i in transcripts) {
                    if (exons.containsKey(i.transcriptId)) {
                        i.exons.addAll(exons[i.transcriptId]!!)
                    }
                }
            }

        }catch (e: IOException) {
            this.logger.error(e.toString())
        }
    }
}