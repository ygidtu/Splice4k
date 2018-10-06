package com.splice4k.index

import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.progressbar.ProgressBar
import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.util.*
import kotlin.system.exitProcess
import com.splice4k.tools.FileValidator


/**
 * 注释文件index的基础类
 * @version 20181006
 * @author Zhang Yiming
 * @since ???
 */


/**
 * @param infile 输入文件
 * @param smrt 是否为SMRT算法运行
 */
class AnnotationIndex(
        private val infile: File,
        private val smrt: Boolean = false,
        private val iso: Boolean = false
) {
    private val logger = Logger.getLogger(AnnotationIndex::class.java)
    private val fileFormat = FileValidator().check(this.infile)
    val data = mutableMapOf<String, MutableList<Exons>>()
    val genes = mutableListOf<Genes>()

    val transcripts: MutableList<Genes> = mutableListOf()

    init {
        this.readFromAnnotation()
    }


    private fun getSourceFromGff(info: List<String>): Map<String, String> {
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


    private fun getSourceFromGtf(info: List<String>): Map<String, String> {
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


    /**
     * 从Gff文件中读取信息
     */
    private fun readFromAnnotation() {


        try {
            val reader = Scanner(this.infile)

            this.logger.info("Reading from ${this.infile}")
            val pb = ProgressBar(message = "Reading exon from Gff")

            val tmpGenes = mutableMapOf<String, Genes>()
            val geneTranscript = mutableMapOf<String, String>()

            val exons: MutableMap<String, MutableList<Exons>> = mutableMapOf()
            var exonNumber = 1

            while (reader.hasNext()) {
                val line = reader.nextLine()

                if ( line.startsWith("#") ) {
                    continue
                }

                val lines = line.replace("\n", "").split("\\s+".toRegex())
                pb.step()

                val sources = when(this.fileFormat) {
                    "gtf" ->  this.getSourceFromGtf(lines.subList(8, lines.size))

                    "gff" -> this.getSourceFromGff(lines.subList(8, lines.size))

                    else -> {
                        logger.info("Please check reference file ${this.infile} format, it should be gtf|gff")
                        exitProcess(2)
                    }
                }


                if ( this.iso && lines[2].matches(".*gene.*".toRegex(RegexOption.IGNORE_CASE)) ) {
                    val gene = Genes(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            information = sources
                    )

                    tmpGenes[gene.geneId] = gene

                }


                if (lines[2].matches(".*exon.*".toRegex(RegexOption.IGNORE_CASE))) {
                    var exonId = sources["exon_id"] ?: sources["ID"]
                    if ( exonId == null ) {
                        exonId = "${sources["Parent"] ?: sources["transcript_id"]}.$exonNumber"
                        exonNumber ++
                    }

                    val tmp = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            exonId = exonId
                    )

                    tmp.source["transcript"]!!.add(sources["Parent"] ?: sources["transcript_id"] ?: sources["ID"]!! )
                    tmp.source["gene"]!!.add(geneTranscript[sources["Parent"]] ?: sources["gene_id"] ?: sources["GeneID"]!!)

                    val tmpExons = mutableListOf(tmp)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmpExons)
                    } else {
                        this.data[key] = tmpExons
                    }

                    if ( this.smrt ) {
                        val tmpExon = mutableListOf(tmp)

                        if (exons.containsKey(sources["Parent"] ?: sources["transcript_id"]!! )) {
                            tmpExon.addAll(exons[sources["Parent"] ?: sources["transcript_id"]]!!)
                        }
                        exons[sources["Parent"] ?: sources["transcript_id"]!!] = tmpExon
                    }

                    if ( this.iso ) {
                        tmpGenes[sources["gene_id"] ?: sources["GeneID"] ?: geneTranscript[sources["Parent"]!!]  ]!!.exons.add(tmp)
                    }
                } else if (  // 两种标准，严防不标准的gff文件
                        ("transcript_id" in sources.keys && "exon_id" !in sources.keys ) ||
                        (lines[2].matches(".*(rna|transcript).*".toRegex(RegexOption.IGNORE_CASE)) &&
                        !lines[2].matches(".*gene.*".toRegex(RegexOption.IGNORE_CASE)))
                ) {

                    geneTranscript[sources["ID"] ?: sources["transcript_id"]!! ] = sources["Parent"] ?: sources["gene_id"] ?: sources["GeneID"]!!

                    if ( this.smrt ) {
                        this.transcripts.add(Genes(
                                chromosome = lines[0],
                                start = lines[3].toInt(),
                                end = lines[4].toInt(),
                                strand = lines[6].toCharArray()[0],
                                information = sources
                        ))
                    }

                    exonNumber = 1
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

            this.genes.addAll(tmpGenes.values)

        }catch (e: IOException) {
            this.logger.error(e.toString())
        }
    }
}